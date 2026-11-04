// Tests du moteur du régulateur PID : node tests/test.js
const assert = require('assert');
const P = require('../pid.js');
let n = 0;
const close = (a, b, tol, msg) => { n++; assert.ok(Math.abs(a - b) <= tol * (1 + Math.abs(b)), (msg || '') + ' : ' + a + ' ≠ ' + b); };
const last = (a) => a[a.length - 1];

// 1. Procédés en boucle ouverte, comparés aux solutions exactes
{
  const p = P.PLANTS.motor.p;
  const s = P.simulate('motor', {}, { gains: {}, openLoop: 12, tEnd: 4, stepAt: 0 });
  close(last(s.y), (12 * p.Kt) / (p.R * p.b + p.Kt * p.Ke), 1e-6, 'moteur : vitesse en régime établi');
  // couveuse : premier ordre avec retard, T(t) = Tamb + K·u·(1 − e^(−(t−θ)/τ))
  const q = P.PLANTS.incubator.p;
  const c = P.simulate('incubator', {}, { gains: {}, openLoop: 0.5, tEnd: 1200, stepAt: 0 });
  [100, 300, 600, 1200].forEach((tt) => {
    const i = c.t.findIndex((v) => v >= tt), want = q.Tamb + q.K * 0.5 * (1 - Math.exp(-Math.max(0, c.t[i] - q.theta) / q.tau));
    close(c.y[i], want, 2e-3, 'couveuse à t = ' + tt + ' s');
  });
  close(c.y[c.t.findIndex((v) => v >= 39)], q.Tamb, 1e-12, 'retard pur : rien ne bouge avant θ');
}

// 2. Fonctions de transfert : gain statique et réponse en fréquence
{
  const p = P.PLANTS.motor.p;
  close(P.cabs(P.plantFreq('motor', {}, 1e-6)), p.Kt / (p.R * p.b + p.Kt * p.Ke), 1e-6, 'gain statique du moteur');
  const q = P.PLANTS.incubator.p, w = 0.01, G = P.plantFreq('incubator', {}, w);
  close(P.cabs(G), q.K / Math.hypot(1, w * q.tau), 1e-12, '|G| couveuse');
  const ph = Math.atan2(G.im, G.re), want = -Math.atan(w * q.tau) - w * q.theta;
  close(Math.atan2(Math.sin(ph - want), Math.cos(ph - want)), 0, 1e-12, 'phase couveuse avec retard');
  // l'axe contient un intégrateur : phase qui part de −90°
  const c = P.bode('axis', {}, null, 1e-3, 1e4, 300);
  close(c[0].phG, -90, 0.01, 'axe : phase basse fréquence −90°');
}

// 3. Gain et période critiques (Ziegler-Nichols) : solution exacte pour la couveuse
{
  const q = P.PLANTS.incubator.p, u = P.ultimate('incubator');
  // −atan(ωτ) − ωθ = −π, puis Ku = √(1 + (ωτ)²)/K
  let lo = 1e-4, hi = 1;
  for (let i = 0; i < 200; i++) { const m = (lo + hi) / 2; if (Math.atan(m * q.tau) + m * q.theta < Math.PI) lo = m; else hi = m; }
  close(u.w, lo, 1e-6, 'pulsation critique'); close(u.Ku, Math.hypot(1, lo * q.tau) / q.K, 1e-6, 'gain critique'); close(u.Tu, (2 * Math.PI) / lo, 1e-6, 'période critique');
  n++; assert.strictEqual(P.ultimate('motor'), null, 'moteur (2e ordre sans retard) : pas de gain critique');
  // au gain critique, la boucle fermée oscille sans s'amortir (marge de gain nulle)
  const g = { Kp: u.Ku, Ki: 0, Kd: 0 };
  close(P.margins(P.bode('incubator', {}, g, 1e-5, 1, 3000)).gm, 0, 0.01, 'marge de gain nulle au gain critique');
  // règles de réglage
  const zn = P.tuning('zn', 10, 2);
  close(zn.Kp, 6, 1e-14, 'ZN Kp'); close(zn.Ki, 6, 1e-14, 'ZN Ki = Kp/(Tu/2)'); close(zn.Kd, 1.5, 1e-14, 'ZN Kd = Kp·Tu/8');
}

// 4. Marges de stabilité sur un cas connu : L(s) = K / (s·(s+1)) a ωc = √((√(1+4K²) − 1)/2)
{
  // axe en P seul : on vérifie la marge de phase sur la courbe calculée, par rapport à la formule de la phase
  const g = { Kp: 12, Ki: 0, Kd: 0 }, c = P.bode('axis', {}, g, 0.01, 1e4, 4000), m = P.margins(c);
  const L = P.cabs(P.plantFreq('axis', {}, m.wc)) * 12;
  close(L, 1, 1e-3, '|L(jωc)| = 1');
  const G = P.plantFreq('axis', {}, m.wc), ph = (Math.atan2(G.im, G.re) * 180) / Math.PI;
  close(m.pm, 180 + (ph > 0 ? ph - 360 : ph), 0.01, 'marge de phase = 180° + arg L(jωc)');
}

// 5. Boucle fermée : erreur statique, action intégrale, perturbation
{
  // couveuse en P seul : erreur statique = (r − Tamb) / (1 + K·Kp) (sans saturation)
  const q = P.PLANTS.incubator.p, Kp = 0.1, r = 30;
  const s = P.simulate('incubator', {}, { gains: { Kp, Ki: 0, Kd: 0 }, setpoint: r, tEnd: 20000, stepAt: 0 });
  close(r - last(s.y), (r - q.Tamb) / (1 + q.K * Kp), 1e-3, 'erreur statique en P seul');
  // avec l'intégrale, l'erreur s'annule, y compris après la perturbation
  const g = Object.assign(P.imcPI(q.K, q.tau, q.theta, 300), { antiWindup: true });
  const t = P.simulate('incubator', {}, { gains: g, tEnd: 9000, dist: true });
  close(last(t.y), P.PLANTS.incubator.setpoint, 1e-3, 'PI : pas d’erreur statique');
  // moteur PI : rejet de la charge
  const m = P.simulate('motor', {}, { gains: { Kp: 0.05, Ki: 2, Kd: 0, antiWindup: true }, dist: true, tEnd: 3 });
  close(last(m.y), 200, 1e-3, 'moteur : charge rejetée');
}

// 6. Saturation et anti-emballement
{
  const sOn = P.simulate('incubator', {}, { gains: { Kp: 0.9, Ki: 0.02, Kd: 0, antiWindup: true } });
  const sOff = P.simulate('incubator', {}, { gains: { Kp: 0.9, Ki: 0.02, Kd: 0, antiWindup: false } });
  n++; assert.ok(sOn.u.every((u) => u >= 0 && u <= 1), 'commande toujours dans [0, 1]');
  const mOn = P.metrics(sOn, 180, 22, 36.5, 3000), mOff = P.metrics(sOff, 180, 22, 36.5, 3000);
  n++; assert.ok(mOff.overshoot > mOn.overshoot + 1, 'sans anti-emballement, le dépassement est plus grand (' + mOff.overshoot.toFixed(1) + ' % contre ' + mOn.overshoot.toFixed(1) + ' %)');
}

// 7. Indicateurs sur une réponse connue (premier ordre : montée 10-90 % = τ·ln 9, réponse à 2 % = τ·ln 50)
{
  const tau = 2, dt = 0.001, t = [], y = [], u = [];
  for (let k = 0; k <= 20000; k++) { t.push(k * dt); y.push(1 - Math.exp(-(k * dt) / tau)); u.push(1); }
  const m = P.metrics({ t, y, u }, 0, 0, 1);
  close(m.rise, tau * Math.log(9), 1e-3, 'temps de montée');
  close(m.settling, tau * Math.log(50), 1e-3, 'temps de réponse à 2 %');
  close(m.overshoot, 0, 1e-12, 'pas de dépassement');
  // second ordre : dépassement = exp(−ζπ/√(1−ζ²))
  const z = 0.3, wn = 5, wd = wn * Math.sqrt(1 - z * z), y2 = t.map((tt) => 1 - Math.exp(-z * wn * tt) * (Math.cos(wd * tt) + (z / Math.sqrt(1 - z * z)) * Math.sin(wd * tt)));
  close(P.metrics({ t, y: y2, u }, 0, 0, 1).overshoot, 100 * Math.exp((-z * Math.PI) / Math.sqrt(1 - z * z)), 1e-3, 'dépassement d’un second ordre');
}

// 8. PID numérique : dérivée sur la mesure (pas de pic de commande au changement de consigne)
{
  const pid = P.makePID({ Kp: 1, Ki: 0, Kd: 1, N: 10 }, 0.01, -1e9, 1e9);
  pid(0, 0);
  const u = pid(1, 0).u; // la consigne saute, la mesure ne bouge pas
  close(u, 1, 1e-12, 'pas de coup de dérivée sur un saut de consigne');
  const v = pid(1, 0.1).u; // la mesure monte : la dérivée freine
  n++; assert.ok(v < 0.9, 'la dérivée s’oppose à la montée de la mesure');
}

console.log('OK : ' + n + ' vérifications passent');
