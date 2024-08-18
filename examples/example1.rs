// mod dubins3d;
use dubins3d::{State, DubinsManeuver3D};
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
    let qi = State{x: 0.0, y: 0.0, z: 0.0, yaw: 0.0, pitch: 0.0};
    let qf = State{x: 100.0, y: 100.0, z: 100.0, yaw: 0.0, pitch: 0.0};
    let rhomin = 10.0;
    let pitchlims = (PI * -15.0 / 180.0, PI * 20.0 / 180.0);

    let dubins = DubinsManeuver3D::new(qi, qf, rhomin, pitchlims);
    let samples = dubins.compute_sampling(500);
    println!("{}", samples.len());

    println!("{} {} {}", samples[0].x - qi.x, samples[0].y - qi.y, samples[0].z - qi.z);
    println!("{} {} {}", samples[samples.len()-1].x - qf.x, samples[samples.len()-1].y - qf.y, samples[samples.len()-1].z - qf.z);
    
    let mut wtr = Writer::from_path("path.csv")?;
    for state in samples.iter() {
        wtr.write_record(&[state.x.to_string(), state.y.to_string(), state.z.to_string()])?;
    }
    wtr.flush()?;
    Ok(())
}
