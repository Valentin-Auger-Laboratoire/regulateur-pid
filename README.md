# Régulateur PID : régler, tester, vérifier la stabilité

🇫🇷 [Français](#français) · 🇬🇧 [English](#english)

**Essayer l'outil :** https://valentin-auger-laboratoire.github.io/portfolio/outils/regulateur-pid/

---

## Français

### Ce que fait l'outil
Simulateur de régulation PID, en HTML / CSS / JavaScript pur (aucune bibliothèque), sur trois systèmes :

| Système | Mesure | Commande | Perturbation |
|---|---|---|---|
| Moteur à courant continu | vitesse (rad/s) | tension, limitée à ±24 V | couple de charge |
| Couveuse de néonatologie (premier ordre avec retard) | température (°C) | puissance de chauffe, de 0 à 100 % | porte ouverte 2 minutes |
| Axe de robot (moteur et réducteur) | angle (°) | tension, limitée à ±24 V | charge sur l'axe |

**Fonctions**
- Gains Kp, Ki et Kd réglables, en échelle logarithmique. Chaque action peut être activée ou désactivée.
- Réponse à un échelon de consigne, puis à la perturbation. Animation du système (disque du moteur, couveuse et thermomètre, bras).
- Indicateurs :
  - dépassement, temps de montée (10 à 90 %), temps de réponse à 2 % ;
  - erreur à la fin, intégrale de l'erreur absolue ;
  - marges de phase et de gain.
- Commande et ses limites, et les trois actions P, I et D séparées.
- **Diagrammes de Bode** du procédé seul et de la boucle ouverte, avec la pulsation de coupure.
- **Réglages automatiques** :
  - Ziegler-Nichols PID et PI, Tyreus-Luyben et la variante « peu de dépassement », à partir du gain critique Ku et de la période critique Tu ;
  - IMC (lambda) pour la couveuse.
  - Quand la méthode ne s'applique pas (moteur : la phase n'atteint jamais −180°), l'outil le dit.
- Anti-emballement de l'intégrale, bruit de mesure, possibilité de figer une courbe pour comparer deux réglages.
- Alertes : instabilité, marge de phase faible, dépassement dangereux, saturation, erreur statique.

### Comment c'est calculé
- Procédés écrits sous forme d'équations d'état et intégrés par **Runge-Kutta 4** entre deux instants d'échantillonnage. Le régulateur est numérique et maintient sa commande constante entre deux calculs.
- PID sous forme parallèle, avec la **dérivée sur la mesure** (pas de coup de commande quand la consigne change) et filtrée (N = 10).
- Anti-emballement par arrêt de l'intégration quand la commande sature.
- Réponse en fréquence calculée sur les fonctions de transfert, avec le retard pur exact e^(−θs). Phase déroulée et calée sur le nombre d'intégrateurs.

### Fichiers
| Fichier | Rôle |
|---|---|
| `pid.js` | Moteur : procédés, PID, simulation, indicateurs, Bode, marges, Ziegler-Nichols, IMC (testable avec Node) |
| `index.html` | Interface, animation, graphiques, bilingue FR / EN |
| `tests/test.js` | Tests du moteur |

### Lancer les tests
```bash
node tests/test.js
```
- Réponses en boucle ouverte comparées aux solutions exactes (moteur, premier ordre avec retard).
- |G| et phase, avec le retard pur, comparés aux formules.
- **Gain et période critiques de la couveuse comparés à la résolution exacte** de −atan(ωτ) − ωθ = −π, et marge de gain nulle au gain critique.
- Erreur statique en P seul égale à (r − T₀)/(1 + K·Kp). Erreur nulle avec l'intégrale, perturbation rejetée.
- L'anti-emballement réduit le dépassement.
- Indicateurs vérifiés sur un premier ordre (τ·ln 9, τ·ln 50) et sur le dépassement d'un second ordre.
- Pas de coup de dérivée sur un saut de consigne.

### Limites
- Paramètres des procédés plausibles, pas ceux d'un matériel précis.
- Pas de frottement sec ni de jeu dans le réducteur de l'axe.
- Modèle de couveuse simplifié (une seule constante de temps, pas d'humidité).
- Les marges sont calculées sur le modèle linéaire, sans la saturation.

### Utiliser en local
Ouvrir `index.html` dans un navigateur.

---

## English

### What the tool does
PID control simulator in plain HTML / CSS / JavaScript, on three systems:
- DC motor speed;
- neonatal incubator temperature (first order with dead time);
- robot axis position.

- Adjustable Kp, Ki, Kd; step response and disturbance (load torque, incubator door opened); animated system.
- Performance figures: overshoot, rise time, settling time, final error, IAE, phase and gain margins.
- Command with its limits and the three P, I, D actions.
- **Bode plots** of the plant and of the open loop.
- **Automatic tunings**: Ziegler-Nichols PID/PI, Tyreus-Luyben and "some overshoot" from Ku and Tu, IMC for the incubator. When the method does not apply, the tool says so.
- Integrator anti-windup, measurement noise, frozen curve to compare tunings, warnings.

### Tests
```bash
node tests/test.js
```
- Open-loop responses versus exact solutions.
- Ultimate gain and period versus the exact solution.
- Steady-state error in P only, zero error with integral action.
- Anti-windup effect; performance figures on first- and second-order responses.
- No derivative kick.

### Limits
- Plausible plant parameters.
- No dry friction or backlash.
- Simplified incubator.
- Margins computed on the linear model.

---

← [Portfolio](https://valentin-auger-laboratoire.github.io/portfolio/#projets)
