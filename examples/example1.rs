// mod dubins3d;
use dubins3d::{DubinsManeuver3D};
// mod vertical;
// mod dubins3d;

use core::f64::consts::PI as PI;
// use crate::dubins3d::{DubinsManeuver3D}; 
use std::error::Error;
use csv::Writer;

// fn example() -> Result<(), Box<dyn Error>> {
//     let mut wtr = Writer::from_path("foo.csv")?;
//     wtr.write_record(&["a", "b", "c"])?;
//     wtr.write_record(&["x", "y", "z"])?;
//     wtr.flush()?;
//     Ok(())
// }

fn main() -> Result<(), Box<dyn Error>> {
    println!("Hello, world!");
    let qi = (0.0, 0.0, 0.0, 0.0, 0.0);
    let qf = (100.0, 100.0, 100.0, 0.0, 0.0);
    let rhomin = 10.0;
    let pitchlims = (PI * -15.0 / 180.0, PI * 20.0 / 180.0);

    let dubins = DubinsManeuver3D::new(qi, qf, rhomin, pitchlims);
    let samples = dubins.compute_sampling(500);
    println!("{}", samples.len());

    println!("{} {} {}", samples[0].0 - qi.0, samples[0].1 - qi.1, samples[0].2 - qi.2);
    println!("{} {} {}", samples[samples.len()-1].0 - qf.0, samples[samples.len()-1].1 - qf.1, samples[samples.len()-1].2 - qf.2);
    
    let mut wtr = Writer::from_path("foo.csv")?;
    for state in samples.iter() {
        wtr.write_record(&[state.0.to_string(), state.1.to_string(), state.2.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}
