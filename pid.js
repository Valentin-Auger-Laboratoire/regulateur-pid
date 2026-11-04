/* Régulateur PID : moteur (sans dépendance).
   Procédés (moteur à courant continu, couveuse, axe de position), PID numérique (dérivée filtrée sur la mesure,
   anti-emballement de l'intégrale, saturation), simulation temporelle, indicateurs de performance,
   réponse en fréquence (Bode), marges de stabilité et réglages de Ziegler-Nichols.
   Utilisable dans le navigateur (window.PID) et avec Node (module.exports). */
(function (root, factory) {
  const M = factory();
  if (typeof module === 'object' && module.exports) module.exports = M;
  else root.PID = M;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---------- Nombres complexes et polynômes ----------
  const cx = (re, im = 0) => ({ re, im });
  const cmul = (a, b) => cx(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
  const cdiv = (a, b) => { const d = b.re * b.re + b.im * b.im; return cx((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d); };
  const cabs = (a) => Math.hypot(a.re, a.im);
  // polynôme en s, coefficients du degré le plus haut au plus bas
  function polyAt(p, s) { return p.reduce((acc, c) => cx(acc.re * s.re - acc.im * s.im + c, acc.re * s.im + acc.im * s.re), cx(0)); }
  const pmul = (a, b) => { const r = new Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => (r[i + j] += x * y))); return r; };
  const padd = (a, b) => { const n = Math.max(a.length, b.length), r = new Array(n).fill(0); a.forEach((x, i) => (r[n - a.length + i] += x)); b.forEach((x, i) => (r[n - b.length + i] += x)); return r; };

  // ---------- Procédés ----------
  // Chaque procédé : paramètres, équations d'état (pour la simulation non linéaire, avec saturation de la commande)
  // et fonction de transfert linéaire G(s) = num/den · e^(−θs) (pour Bode et les réglages)
  const PLANTS = {
    // Moteur à courant continu, sortie = vitesse (rad/s), commande = tension (V)
    motor: {
      p: { J: 2e-4, b: 1e-4, Kt: 0.05, Ke: 0.05, R: 2, L: 1e-3, umin: -24, umax: 24, load: 0.02 },
      x0: () => [0, 0], // [ω, i]
      deriv(x, u, d, p) { const [w, i] = x; return [(p.Kt * i - p.b * w - d * p.load) / p.J, (u - p.R * i - p.Ke * w) / p.L]; },
      y: (x) => x[0],
      tf(p) { return { num: [p.Kt], den: padd(pmul([p.L, p.R], [p.J, p.b]), [p.Kt * p.Ke]), delay: 0 }; },
      setpoint: 200, tEnd: 1.5, Ts: 0.001, dist: { t: 0.8, label: 'load' }, y0: 0,
    },
    // Couveuse : premier ordre avec retard, sortie = température (°C), commande = puissance de chauffe (0 à 1)
    incubator: {
      p: { K: 16, tau: 600, theta: 40, Tamb: 22, umin: 0, umax: 1, loss: 5 },
      x0: (p) => [p.Tamb],
      // la porte ouverte (d = 1) ajoute une perte équivalente à « loss » °C
      deriv(x, u, d, p) { return [(-(x[0] - p.Tamb) + p.K * u - d * p.loss) / p.tau]; },
      y: (x) => x[0],
      tf(p) { return { num: [p.K], den: [p.tau, 1], delay: p.theta }; },
      setpoint: 36.5, tEnd: 3600, Ts: 1, dist: { t: 1800, dur: 120, label: 'door' }, y0: 22,
    },
    // Axe de position (moteur + réducteur), sortie = angle (°), commande = tension (V)
    axis: {
      p: { J: 2e-4, b: 1e-4, Kt: 0.05, Ke: 0.05, R: 2, L: 1e-3, N: 50, Jl: 0.2, umin: -24, umax: 24, load: 1.5 },
      x0: () => [0, 0, 0], // [angle de sortie (rad), vitesse moteur (rad/s), i]
      deriv(x, u, d, p) {
        const [th, w, i] = x, Jeq = p.J + p.Jl / (p.N * p.N);
        return [w / p.N, (p.Kt * i - p.b * w - (d * p.load) / p.N) / Jeq, (u - p.R * i - p.Ke * w) / p.L];
      },
      y: (x) => (x[0] * 180) / Math.PI,
      tf(p) {
        const Jeq = p.J + p.Jl / (p.N * p.N);
        // θ(s)/U(s) en degrés : (180/π)·Kt / (N·s·((Ls+R)(Jeq·s+b) + Kt·Ke))
        return { num: [((180 / Math.PI) * p.Kt) / p.N], den: pmul([1, 0], padd(pmul([p.L, p.R], [Jeq, p.b]), [p.Kt * p.Ke])), delay: 0 };
      },
      setpoint: 90, tEnd: 3, Ts: 0.001, dist: { t: 1.6, label: 'load' }, y0: 0,
    },
  };

  // ---------- Contrôleur PID numérique ----------
  // forme parallèle : u = Kp·e + Ki·∫e + Kd·dy_f/dt (dérivée sur la mesure, filtrée à Kd/(Kp·N))
  function makePID(g, Ts, umin, umax) {
    let I = 0, yPrev = null, D = 0;
    const N = g.N || 10, Tf = g.Kp > 0 && g.Kd > 0 ? g.Kd / (g.Kp * N) : 0;
    return function step(r, y) {
      const e = r - y;
      if (yPrev === null) yPrev = y;
      // dérivée filtrée (filtre du premier ordre, discrétisé par Euler implicite)
      if (g.Kd > 0) D = (Tf * D - g.Kd * (y - yPrev)) / (Tf + Ts);
      yPrev = y;
      let u = g.Kp * e + I + D;
      const sat = Math.min(umax, Math.max(umin, u));
      // anti-emballement : on n'intègre pas si la commande est saturée et que l'erreur pousse dans le même sens
      if (!(g.antiWindup && sat !== u && Math.sign(e) === Math.sign(u - sat))) I += g.Ki * e * Ts;
      return { u: sat, uRaw: u };
    };
  }

  // ---------- Simulation en boucle fermée ----------
  // opts : { gains, setpoint, tEnd, Ts, dist (bool), noise (écart-type de mesure), sub (sous-pas RK4), openLoop: u fixe }
  function simulate(key, params, opts) {
    const P = PLANTS[key], p = Object.assign({}, P.p, params || {});
    const Ts = opts.Ts || P.Ts, tEnd = opts.tEnd || P.tEnd, r0 = opts.setpoint === undefined ? P.setpoint : opts.setpoint;
    const sub = opts.sub || 4, h = Ts / sub;
    const delaySteps = key === 'incubator' ? Math.round(p.theta / Ts) : 0;
    const uBuf = new Array(delaySteps).fill(0);
    const pid = makePID(opts.gains, Ts, p.umin, p.umax);
    let x = P.x0(p), t = 0;
    const out = { t: [], y: [], u: [], r: [], d: [] };
    let seed = 12345;
    const gauss = () => { seed = (seed * 16807) % 2147483647; const u1 = seed / 2147483647; seed = (seed * 16807) % 2147483647; const u2 = seed / 2147483647; return Math.sqrt(-2 * Math.log(u1 + 1e-12)) * Math.cos(2 * Math.PI * u2); };
    const stepAt = opts.stepAt === undefined ? tEnd * 0.05 : opts.stepAt;
    const n = Math.round(tEnd / Ts);
    for (let k = 0; k <= n; k++) {
      const r = t >= stepAt ? r0 : P.y0;
      const dOn = opts.dist && P.dist && t >= P.dist.t && (!P.dist.dur || t < P.dist.t + P.dist.dur) ? 1 : 0;
      const ym = P.y(x) + (opts.noise ? opts.noise * gauss() : 0);
      const u = opts.openLoop !== undefined ? opts.openLoop : pid(r, ym).u;
      out.t.push(t); out.y.push(P.y(x)); out.u.push(u); out.r.push(r); out.d.push(dOn);
      // retard pur sur la commande (couveuse)
      let ua = u;
      if (delaySteps) { uBuf.push(u); ua = uBuf.shift(); }
      for (let s = 0; s < sub; s++) x = rk4(P.deriv, x, ua, dOn, p, h);
      t += Ts;
    }
    return out;
  }
  function rk4(f, x, u, d, p, h) {
    const k1 = f(x, u, d, p), a = x.map((v, i) => v + (h / 2) * k1[i]);
    const k2 = f(a, u, d, p), b = x.map((v, i) => v + (h / 2) * k2[i]);
    const k3 = f(b, u, d, p), c = x.map((v, i) => v + h * k3[i]);
    const k4 = f(c, u, d, p);
    return x.map((v, i) => v + (h / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
  }

  // ---------- Indicateurs de la réponse à un échelon ----------
  function metrics(sim, stepAt, y0, r, tUntil) {
    const { t, y, u } = sim, A = r - y0, dir = Math.sign(A) || 1;
    let i0 = t.findIndex((v) => v >= stepAt), i1 = tUntil ? t.findIndex((v) => v >= tUntil) : t.length;
    if (i1 < 0) i1 = t.length;
    const yy = y.slice(i0, i1), tt = t.slice(i0, i1);
    const rel = yy.map((v) => (v - y0) / A);
    const peak = Math.max(...rel);
    const overshoot = Math.max(0, (peak - 1) * 100);
    const t10 = tt[rel.findIndex((v) => v >= 0.1)], t90 = tt[rel.findIndex((v) => v >= 0.9)];
    let ts = null;
    for (let k = rel.length - 1; k >= 0; k--) if (Math.abs(rel[k] - 1) > 0.02) { ts = k + 1 < tt.length ? tt[k + 1] - stepAt : null; break; }
    if (ts === null && Math.abs(rel[0] - 1) <= 0.02) ts = 0;
    const ess = r - yy[yy.length - 1];
    let iae = 0;
    for (let k = 1; k < yy.length; k++) iae += Math.abs(r - yy[k]) * (tt[k] - tt[k - 1]);
    const umax = Math.max(...u.slice(i0, i1).map(Math.abs));
    return { overshoot, rise: t10 !== undefined && t90 !== undefined ? t90 - t10 : null, settling: ts, ess, iae, umax, peak: y0 + peak * A, dir };
  }

  // ---------- Réponse en fréquence ----------
  // G(jω) du procédé, C(jω) du PID (avec le filtre de dérivée), L = C·G
  function plantFreq(key, params, w) {
    const p = Object.assign({}, PLANTS[key].p, params || {}), tf = PLANTS[key].tf(p);
    const s = cx(0, w);
    let G = cdiv(polyAt(tf.num, s), polyAt(tf.den, s));
    if (tf.delay) G = cmul(G, cx(Math.cos(-w * tf.delay), Math.sin(-w * tf.delay)));
    return G;
  }
  function pidFreq(g, w) {
    const N = g.N || 10, Tf = g.Kp > 0 && g.Kd > 0 ? g.Kd / (g.Kp * N) : 0, s = cx(0, w);
    let C = cx(g.Kp);
    if (g.Ki) C = cx(C.re + cdiv(cx(g.Ki), s).re, C.im + cdiv(cx(g.Ki), s).im);
    if (g.Kd) { const D = cdiv(cmul(cx(g.Kd), s), cx(1, w * Tf)); C = cx(C.re + D.re, C.im + D.im); }
    return C;
  }
  // phase déroulée sur une grille de pulsations croissantes
  function bode(key, params, g, w0, w1, n = 400) {
    const out = [];
    let prev = null, off = 0, prevP = null, offP = 0;
    for (let k = 0; k <= n; k++) {
      const w = w0 * Math.pow(w1 / w0, k / n);
      const G = plantFreq(key, params, w), L = g ? cmul(pidFreq(g, w), G) : G;
      let ph = Math.atan2(L.im, L.re), phG = Math.atan2(G.im, G.re);
      if (prev !== null) { while (ph + off - prev > Math.PI) off -= 2 * Math.PI; while (ph + off - prev < -Math.PI) off += 2 * Math.PI; }
      if (prevP !== null) { while (phG + offP - prevP > Math.PI) offP -= 2 * Math.PI; while (phG + offP - prevP < -Math.PI) offP += 2 * Math.PI; }
      prev = ph + off; prevP = phG + offP;
      out.push({ w, mag: 20 * Math.log10(cabs(L)), ph: (prev * 180) / Math.PI, magG: 20 * Math.log10(cabs(G)), phG: (prevP * 180) / Math.PI });
    }
    // la phase n'est connue qu'à 360° près : on cale le premier point sur la valeur attendue aux basses
    // fréquences, −90° par intégrateur (pôles en 0 du procédé, plus l'action intégrale du PID)
    const tf = PLANTS[key].tf(Object.assign({}, PLANTS[key].p, params || {}));
    let k0 = 0; for (let i = tf.den.length - 1; i >= 0 && tf.den[i] === 0; i--) k0++;
    const fix = (field, nInt) => {
      const want = -90 * nInt, d = Math.round((out[0][field] - want) / 360) * 360;
      if (d) out.forEach((q) => (q[field] -= d));
    };
    fix('phG', k0);
    fix('ph', k0 + (g && g.Ki > 0 ? 1 : 0));
    return out;
  }
  // Marges : de phase (à |L| = 1) et de gain (à phase = −180°)
  function margins(curve) {
    let pm = null, wc = null, gm = null, w180 = null;
    for (let k = 1; k < curve.length; k++) {
      const a = curve[k - 1], b = curve[k];
      if (pm === null && a.mag >= 0 && b.mag < 0) { const u = a.mag / (a.mag - b.mag); wc = a.w * Math.pow(b.w / a.w, u); const ph = a.ph + (b.ph - a.ph) * u; pm = 180 + ph; }
      if (gm === null && a.ph > -180 && b.ph <= -180) { const u = (a.ph + 180) / (a.ph - b.ph); w180 = a.w * Math.pow(b.w / a.w, u); const mag = a.mag + (b.mag - a.mag) * u; gm = -mag; }
    }
    return { pm, wc, gm, w180 };
  }
  // Gain critique et période critique du procédé seul (pour Ziegler-Nichols) : phase de G = −180°
  function ultimate(key, params) {
    const tf = PLANTS[key].tf(Object.assign({}, PLANTS[key].p, params || {}));
    const c = bode(key, params, null, 1e-5, 1e5, 3000);
    for (let k = 1; k < c.length; k++) {
      const a = c[k - 1], b = c[k];
      if (a.phG > -180 && b.phG <= -180) {
        // affinage par dichotomie sur la phase
        let lo = a.w, hi = b.w;
        const phaseAt = (w) => { const cc = bode(key, params, null, w * 1e-6, w, 2000); return cc[cc.length - 1].phG; };
        for (let i = 0; i < 60; i++) { const m = Math.sqrt(lo * hi); if (phaseAt(m) > -180) lo = m; else hi = m; }
        const w = Math.sqrt(lo * hi), G = plantFreq(key, params, w);
        return { Ku: 1 / cabs(G), Tu: (2 * Math.PI) / w, w };
      }
    }
    return null; // la phase n'atteint jamais −180° : Ziegler-Nichols (méthode du pompage) ne s'applique pas
  }
  // Règles de réglage à partir de Ku et Tu
  function tuning(rule, Ku, Tu) {
    let Kp, Ti, Td;
    if (rule === 'zn') { Kp = 0.6 * Ku; Ti = Tu / 2; Td = Tu / 8; }
    else if (rule === 'zn-pi') { Kp = 0.45 * Ku; Ti = Tu / 1.2; Td = 0; }
    else if (rule === 'tl') { Kp = 0.45 * Ku; Ti = 2.2 * Tu; Td = Tu / 6.3; } // Tyreus-Luyben, plus doux
    else if (rule === 'nos') { Kp = 0.2 * Ku; Ti = Tu / 2; Td = Tu / 3; } // « sans dépassement »
    return { Kp, Ki: Ti ? Kp / Ti : 0, Kd: Kp * Td };
  }
  // Réglage IMC (lambda) pour un premier ordre avec retard K, τ, θ : régulateur PI
  function imcPI(K, tau, theta, lambda) { const Kp = tau / (K * (lambda + theta)); return { Kp, Ki: Kp / tau, Kd: 0 }; }

  return { PLANTS, makePID, simulate, metrics, plantFreq, pidFreq, bode, margins, ultimate, tuning, imcPI, polyAt, pmul, padd, cx, cabs };
});
