use rand::Rng;

/// 2D ferromagnetic Ising model on an L×L square lattice
/// with periodic boundary conditions. J = k_B = 1.
pub struct IsingModel {
    pub l: usize,
    pub spins: Vec<i8>, // +1 or -1
    // Precomputed Metropolis acceptance probabilities.
    // ΔE = 2 · s_i · Σ_nn s_j; since each spin is ±1 and there are 4
    // neighbours, s_i · Σ_nn ∈ {−4,−2,0,2,4}.
    // Index k = (s_i · Σ_nn + 4) / 2 ∈ {0,1,2,3,4}
    // exp_table[k] = min(1, exp(−ΔE / T))
    exp_table: [f64; 5],
}

impl IsingModel {
    /// All spins up – good starting point for T < T_c.
    pub fn new_ordered(l: usize, t: f64) -> Self {
        let spins = vec![1i8; l * l];
        Self::from_spins(l, t, spins)
    }

    /// Random initial configuration.
    pub fn new_random<R: Rng>(l: usize, t: f64, rng: &mut R) -> Self {
        let spins: Vec<i8> = (0..l * l)
            .map(|_| if rng.gen_bool(0.5) { 1 } else { -1 })
            .collect();
        Self::from_spins(l, t, spins)
    }

    fn from_spins(l: usize, t: f64, spins: Vec<i8>) -> Self {
        let mut model = IsingModel { l, spins, exp_table: [0.0; 5] };
        model.set_temperature(t);
        model
    }

    /// Update the Metropolis table when scanning temperatures.
    pub fn set_temperature(&mut self, t: f64) {
        self.exp_table = [
            1.0,               // s·Σ = −4, ΔE = −8: always accept
            1.0,               // s·Σ = −2, ΔE = −4: always accept
            1.0,               // s·Σ =  0, ΔE =  0: always accept
            (-4.0 / t).exp(), // s·Σ =  2, ΔE =  4
            (-8.0 / t).exp(), // s·Σ =  4, ΔE =  8
        ];
    }

    /// Sum of the four nearest-neighbour spins at lattice site i.
    #[inline]
    fn neighbor_sum(&self, i: usize) -> i32 {
        let l = self.l;
        let r = i / l;
        let c = i % l;
        self.spins[((r + l - 1) % l) * l + c] as i32  // up
            + self.spins[((r + 1) % l) * l + c] as i32 // down
            + self.spins[r * l + (c + l - 1) % l] as i32 // left
            + self.spins[r * l + (c + 1) % l] as i32     // right
    }

    /// One Monte Carlo sweep: L² random single-spin flip attempts.
    pub fn sweep<R: Rng>(&mut self, rng: &mut R) {
        let n = self.l * self.l;
        for _ in 0..n {
            let i = rng.gen_range(0..n);
            let s = self.spins[i] as i32;
            let sum_nn = self.neighbor_sum(i);
            let idx = ((s * sum_nn + 4) / 2) as usize;
            if self.exp_table[idx] >= 1.0 || rng.gen::<f64>() < self.exp_table[idx] {
                self.spins[i] = -self.spins[i];
            }
        }
    }

    /// Total energy E = −J Σ_{⟨i,j⟩} s_i s_j.
    /// Counts right + down bonds at every site (each bond counted once for L ≥ 3).
    /// NOTE: for L = 2 with periodic BC every bond is double-counted; the value
    /// is still self-consistent with delta_energy, but is 2× the physical energy.
    pub fn energy(&self) -> f64 {
        let l = self.l;
        let mut e = 0i32;
        for r in 0..l {
            for c in 0..l {
                let s = self.spins[r * l + c] as i32;
                e -= s * self.spins[r * l + (c + 1) % l] as i32; // right
                e -= s * self.spins[((r + 1) % l) * l + c] as i32; // down
            }
        }
        e as f64
    }

    /// Energy change when flipping spin at site i (without actually flipping).
    pub fn delta_energy(&self, i: usize) -> f64 {
        2.0 * self.spins[i] as f64 * self.neighbor_sum(i) as f64
    }

    /// Raw magnetization sum M = Σ s_i.
    pub fn magnetization_sum(&self) -> f64 {
        self.spins.iter().map(|&s| s as f64).sum()
    }
}

// ─── Unit tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{rngs::SmallRng, SeedableRng};

    fn all_up(l: usize) -> IsingModel {
        IsingModel::new_ordered(l, 2.0)
    }

    fn all_down(l: usize) -> IsingModel {
        let mut m = IsingModel::new_ordered(l, 2.0);
        m.spins.iter_mut().for_each(|s| *s = -1);
        m
    }

    fn checkerboard(l: usize) -> IsingModel {
        let mut m = IsingModel::new_ordered(l, 2.0);
        for r in 0..l {
            for c in 0..l {
                if (r + c) % 2 == 1 {
                    m.spins[r * l + c] = -1;
                }
            }
        }
        m
    }

    // ── energy (tested on L = 4 where each bond is counted exactly once) ──

    #[test]
    fn energy_all_up_l4() {
        // 2 × L² bonds (right + down), all s_i s_j = +1 → E = −2L² = −32
        assert_eq!(all_up(4).energy(), -32.0);
    }

    #[test]
    fn energy_all_down_l4() {
        // same bond structure, all s_i s_j = +1 → E = −32
        assert_eq!(all_down(4).energy(), -32.0);
    }

    #[test]
    fn energy_checkerboard_l4() {
        // all bonds are between opposite spins → E = +32
        assert_eq!(checkerboard(4).energy(), 32.0);
    }

    // ── ΔE ──────────────────────────────────────────────────────────────────

    #[test]
    fn delta_energy_all_up_l4() {
        // Flipping any spin costs 4 bonds × 2 = 8
        let m = all_up(4);
        assert_eq!(m.delta_energy(0), 8.0);
        assert_eq!(m.delta_energy(7), 8.0);
        assert_eq!(m.delta_energy(15), 8.0);
    }

    #[test]
    fn delta_energy_checkerboard_l4() {
        // All 4 neighbours are opposite → ΔE = −8
        let m = checkerboard(4);
        assert_eq!(m.delta_energy(0), -8.0);
        assert_eq!(m.delta_energy(1), -8.0);
    }

    #[test]
    fn delta_energy_consistent_with_energy() {
        // delta_energy(i) must equal energy_after − energy_before
        let mut m = all_up(4);
        for i in [0usize, 1, 5, 10, 15] {
            let e_before = m.energy();
            let de = m.delta_energy(i);
            m.spins[i] = -m.spins[i];
            let e_after = m.energy();
            m.spins[i] = -m.spins[i]; // restore
            assert!(
                (e_after - e_before - de).abs() < 1e-10,
                "site {i}: delta_energy = {de}, actual ΔE = {}",
                e_after - e_before
            );
        }
    }

    // ── magnetization ───────────────────────────────────────────────────────

    #[test]
    fn magnetization_all_up_l4() {
        assert_eq!(all_up(4).magnetization_sum(), 16.0);
    }

    #[test]
    fn magnetization_all_down_l4() {
        assert_eq!(all_down(4).magnetization_sum(), -16.0);
    }

    #[test]
    fn magnetization_all_up_l2() {
        assert_eq!(all_up(2).magnetization_sum(), 4.0);
    }

    // ── physical limits ─────────────────────────────────────────────────────

    #[test]
    fn low_temperature_stays_ordered() {
        // At T = 0.5 ≪ T_c ≈ 2.269, the ordered state should not disorder.
        let mut rng = SmallRng::seed_from_u64(42);
        let l = 4;
        let mut m = IsingModel::new_ordered(l, 0.5);
        for _ in 0..10_000 {
            m.sweep(&mut rng);
        }
        let mag = m.magnetization_sum().abs() / (l * l) as f64;
        assert!(mag > 0.99, "low-T |m| = {mag:.4}, expected ≈ 1.0");
    }

    #[test]
    fn high_temperature_mean_magnetization_near_zero() {
        // At T = 50 ≫ T_c, time-averaged ⟨m⟩ should be ≈ 0.
        let mut rng = SmallRng::seed_from_u64(123);
        let l = 4;
        let n = (l * l) as f64;
        let mut model = IsingModel::new_random(l, 50.0, &mut rng);
        for _ in 0..2_000 { // thermalize
            model.sweep(&mut rng);
        }
        let mut sum_m = 0.0f64;
        let steps = 50_000usize;
        for _ in 0..steps {
            model.sweep(&mut rng);
            sum_m += model.magnetization_sum() / n;
        }
        let mean_m = sum_m / steps as f64;
        assert!(
            mean_m.abs() < 0.05,
            "high-T ⟨m⟩ = {mean_m:.4}, expected ≈ 0.0"
        );
    }
}
