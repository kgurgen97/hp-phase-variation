"""Technical-noise (sequencing error + repeat stutter) model and its estimation from pure-allele data."""
import math
from collections import defaultdict

from .common import motif_class


def binom_sf(x, n, p):
    """P(X >= x) for X ~ Binomial(n, p), log-space, exact."""
    if x <= 0:
        return 1.0
    if x > n:
        return 0.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    lp, lq = math.log(p), math.log1p(-p)
    lg = math.lgamma
    terms = []
    for i in range(x, n + 1):
        terms.append(lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * lp + (n - i) * lq)
        if len(terms) > 3 and terms[-1] < terms[0] - 40 and i > n * p:
            break
    mx = max(terms)
    return min(1.0, math.exp(mx) * sum(math.exp(t - mx) for t in terms))


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 1.0
    ph = k / n
    d = 1 + z * z / n
    c = ph + z * z / (2 * n)
    a = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return max(0.0, (c - a) / d), min(1.0, (c + a) / d)


class NoiseModel:
    """rate(k units from dominant): |k|=1 stutter a+b*copies (down) or *up_ratio (up); |k|=2 times k2_ratio;
    larger distances or none-of-the-above: background_rate. All parameters live in the rules table."""

    def __init__(self, params):
        self.p = params

    def stutter1(self, mclass, copies, direction):
        a = self.p["stutter_a_" + mclass]
        b = self.p["stutter_b_" + mclass]
        r = max(0.0, min(0.5, a + b * copies))
        return r if direction < 0 else r * self.p.stutter_up_ratio

    def rate(self, mclass, k, copies):
        """Expected per-read fraction of an allele k motif units away from the dominant allele."""
        if k == 0:
            return 0.0
        d = -1 if k < 0 else 1
        ak = abs(k)
        if ak <= 2:
            r = self.stutter1(mclass, copies, d) * (self.p.stutter_k2_ratio ** (ak - 1))
            return max(r, self.p.background_rate)
        return self.p.background_rate


def _wls(xs, ys, ws):
    sw = sum(ws)
    if sw <= 0:
        return None
    mx = sum(w * x for w, x in zip(ws, xs)) / sw
    my = sum(w * y for w, y in zip(ws, ys)) / sw
    sxx = sum(w * (x - mx) ** 2 for w, x in zip(ws, xs))
    if sxx < 1e-9:
        return my, 0.0
    b = sum(w * (x - mx) * (y - my) for w, x, y in zip(ws, xs, ys)) / sxx
    return my - b * mx, b


def estimate_noise(observations, params, min_obs=3):
    """Estimate stutter parameters from PURE-allele sample x locus observations.

    observations: dicts {mclass, copies (dominant), n (unit-aligned reads), m1, p1 (reads at -1/+1 unit),
    m2, p2 (-2/+2), far (>=3 units), offunit (non unit-aligned reads), n_all (aligned + offunit)}.
    Per motif class: weighted least squares of the -1 unit rate on dominant copies (weights = n), clipped
    to >= 0. up_ratio = sum(p1)/sum(m1); k2_ratio = sum(m2+p2)/sum(m1+p1); background_rate = (far+offunit)/n_all.
    Returns (dict of param_id -> value, report list)."""
    out, report = {}, []
    by = defaultdict(list)
    for o in observations:
        by[o["mclass"]].append(o)
    for mc in ("MONO", "DI", "POLY"):
        obs = [o for o in by.get(mc, []) if o["n"] > 0]
        if len(obs) < min_obs:
            report.append("%s: %d observations < %d, defaults retained" % (mc, len(obs), min_obs))
            continue
        fit = _wls([o["copies"] for o in obs], [o["m1"] / o["n"] for o in obs], [o["n"] for o in obs])
        a, b = fit
        if b < 0:
            b = 0.0
            a = sum(o["m1"] for o in obs) / sum(o["n"] for o in obs)
        out["stutter_a_" + mc] = round(max(a, 0.0), 6)
        out["stutter_b_" + mc] = round(b, 6)
        report.append("%s: n_obs=%d a=%.5f b=%.5f" % (mc, len(obs), a, b))
    sm1 = sum(o["m1"] for o in observations)
    sp1 = sum(o["p1"] for o in observations)
    if sm1 >= 10:
        out["stutter_up_ratio"] = round(sp1 / sm1, 4)
    s1 = sm1 + sp1
    s2 = sum(o["m2"] + o["p2"] for o in observations)
    if s1 >= 10:
        out["stutter_k2_ratio"] = round(s2 / s1, 4)
    na = sum(o.get("n_all", o["n"]) for o in observations)
    if na > 0:
        # floored at 0.001 so that a small pure-allele data set cannot claim a zero error floor
        out["background_rate"] = round(max(0.001, (sum(o["far"] + o["offunit"] for o in observations)) / na), 6)
    return out, report
