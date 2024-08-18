pub(crate) use core::f64::consts::PI as PI;
pub(crate) fn mod2pi(th: f64) -> f64{
    let t = th % (2.0 * PI);
    if t < 0.0 {
        return t + 2.0 * PI;
    }
    return t;
}

pub struct state {
    x: f64,
    y: f64,
    z: f64,
    yaw: f64,
    pitch: f64
}

mod dubins2d;
mod vertical;
mod dubins3d;

pub use crate::dubins3d::{DubinsManeuver3D};