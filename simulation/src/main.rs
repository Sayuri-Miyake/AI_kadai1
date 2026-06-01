mod ising;
mod observables;

use std::io::{self, Write};
use std::time::Instant;

use rand::{rngs::SmallRng, SeedableRng};

use ising::IsingModel;
use observables::Accumulator;

// ─── Simulation parameters ──────────────────────────────────────────────────

const SIZES: &[usize] = &[8, 16, 32, 64];

// Temperature grid: coarse scan + finer grid near T_c for finite-size scaling
const T_COARSE_MIN: f64 = 1.60;
const T_COARSE_MAX: f64 = 2.10;
const T_COARSE_STEP: f64 = 0.10;

const T_FINE_MIN: f64 = 2.10;
const T_FINE_MAX: f64 = 2.45;
const T_FINE_STEP: f64 = 0.02;

const T_HIGH_MIN: f64 = 2.45;
const T_HIGH_MAX: f64 = 3.00;
const T_HIGH_STEP: f64 = 0.10;

// Onsager exact result T_c / J = 2 / ln(1 + √2)
const T_C_ONSAGER: f64 = 2.269_185_314_213_022;

// ─── Sweep counts ───────────────────────────────────────────────────────────

fn sweep_counts(l: usize) -> (usize, usize) {
    // (thermalization MCS, measurement MCS)
    // Autocorrelation time grows as L^z with z ≈ 2 near T_c for local updates.
    match l {
        0..=16 => (1_000, 10_000),
        17..=32 => (3_000, 30_000),
        _ => (8_000, 80_000),
    }
}

// ─── Temperature grid ───────────────────────────────────────────────────────

fn temperature_grid() -> Vec<f64> {
    let mut temps = Vec::new();
    let ranges = [
        (T_COARSE_MIN, T_COARSE_MAX, T_COARSE_STEP),
        (T_FINE_MIN, T_FINE_MAX, T_FINE_STEP),
        (T_HIGH_MIN, T_HIGH_MAX, T_HIGH_STEP),
    ];
    for (tmin, tmax, step) in ranges {
        let mut t = tmin;
        while t <= tmax + 1e-10 {
            temps.push(t);
            t += step;
        }
    }
    temps.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
    temps
}

// ─── Main simulation loop ───────────────────────────────────────────────────

fn main() {
    let mut rng = SmallRng::from_entropy();
    let temps = temperature_grid();
    let stdout = io::stdout();
    let mut out = stdout.lock();

    writeln!(out, "L,T,m,chi,C,U_L").unwrap();
    eprintln!("2D Ising MC  (T_c_Onsager = {T_C_ONSAGER:.6})");
    eprintln!("{:>4}  {:>6}  {:>8}  {:>8}  {:>8}  {:>8}",
        "L", "T", "m", "chi", "C", "U_L");

    for &l in SIZES {
        let n = l * l;
        let (n_therm, n_meas) = sweep_counts(l);
        let t0 = Instant::now();

        // Build an initial ordered lattice, then scan temperatures in order.
        // Reusing the previous equilibrium configuration as the starting point
        // for each new T significantly reduces the required thermalization.
        let mut model = IsingModel::new_ordered(l, temps[0]);

        for &t in &temps {
            model.set_temperature(t);

            for _ in 0..n_therm {
                model.sweep(&mut rng);
            }

            let mut acc = Accumulator::new();
            for _ in 0..n_meas {
                model.sweep(&mut rng);
                let e = model.energy() / n as f64;
                let m = model.magnetization_sum() / n as f64;
                acc.push(e, m);
            }

            let obs = acc.compute(n, t);
            writeln!(
                out,
                "{},{:.4},{:.6},{:.6},{:.6},{:.6}",
                l, t, obs.m, obs.chi, obs.specific_heat, obs.binder
            )
            .unwrap();

            eprintln!(
                "{:>4}  {:>6.3}  {:>8.4}  {:>8.3}  {:>8.3}  {:>8.4}",
                l, t, obs.m, obs.chi, obs.specific_heat, obs.binder
            );
        }

        eprintln!("  → L={l} done in {:.1}s", t0.elapsed().as_secs_f64());
    }
}
