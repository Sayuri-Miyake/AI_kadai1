/// Accumulates raw moments of energy-per-spin and magnetization-per-spin
/// over measurement sweeps, then computes derived physical observables.
pub struct Accumulator {
    n: u64,
    sum_e: f64,
    sum_e2: f64,
    sum_abm: f64, // Σ |m|
    sum_m2: f64,  // Σ m²
    sum_m4: f64,  // Σ m⁴
}

pub struct Observables {
    /// Mean absolute magnetization per spin ⟨|m|⟩
    pub m: f64,
    /// Magnetic susceptibility per spin χ = (N/T)(⟨m²⟩ − ⟨|m|⟩²)
    pub chi: f64,
    /// Specific heat per spin C = (N/T²)(⟨e²⟩ − ⟨e⟩²)
    pub specific_heat: f64,
    /// Binder cumulant U_L = 1 − ⟨m⁴⟩ / (3⟨m²⟩²)
    pub binder: f64,
}

impl Accumulator {
    pub fn new() -> Self {
        Accumulator {
            n: 0,
            sum_e: 0.0,
            sum_e2: 0.0,
            sum_abm: 0.0,
            sum_m2: 0.0,
            sum_m4: 0.0,
        }
    }

    /// Record one sample.
    /// `e` = energy per spin (E_total / N), `m` = magnetization per spin (M / N).
    pub fn push(&mut self, e: f64, m: f64) {
        self.n += 1;
        self.sum_e += e;
        self.sum_e2 += e * e;
        let abm = m.abs();
        let m2 = m * m;
        self.sum_abm += abm;
        self.sum_m2 += m2;
        self.sum_m4 += m2 * m2;
    }

    /// Compute physical observables from accumulated moments.
    /// `n_sites` = L², `t` = temperature.
    pub fn compute(&self, n_sites: usize, t: f64) -> Observables {
        assert!(self.n > 0, "no samples accumulated");
        let n = self.n as f64;
        let vol = n_sites as f64;

        let mean_e = self.sum_e / n;
        let mean_e2 = self.sum_e2 / n;
        let mean_abm = self.sum_abm / n;
        let mean_m2 = self.sum_m2 / n;
        let mean_m4 = self.sum_m4 / n;

        let chi = vol / t * (mean_m2 - mean_abm * mean_abm);
        let specific_heat = vol / (t * t) * (mean_e2 - mean_e * mean_e);
        let binder = if mean_m2 < 1e-30 {
            0.0
        } else {
            1.0 - mean_m4 / (3.0 * mean_m2 * mean_m2)
        };

        Observables { m: mean_abm, chi, specific_heat, binder }
    }
}
