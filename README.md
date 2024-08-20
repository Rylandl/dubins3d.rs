# Dubins3d.rs

This library is a Rust reimplementation of Petr Vana's [Dubins3D.jl](https://github.com/comrob/Dubins3D.jl/tree/master) Julia library. It allows the generation of 3D Dubins paths between two states while bounded by turn radius and pitch angle constraits.

## Example
An example is provided with the repo that generates a route and writes the generated path to `path.csv` file. 

```bash
cargo run --example example1
```
The result can then be visualized with the provided python script. This will save the plot as `path.png`

```bash
python3 scripts/plot_path.py path.csv
```

## Usage
```rust
use core::f64::consts::PI as PI;
use dubins3d::{State, DubinsManeuver3D};

fn main() {
    let qi = State{
        x: 0.0, y: 0.0, z: 0.0, yaw: 0.0, pitch: 0.0
    };
    let qf = State{
        x: 100.0, y: 100.0, z: 100.0, yaw: 0.0, pitch: 0.0
    };
    let min_turn_radius = 10.0;
    let pitch_lims = (PI * -15.0 / 180.0, PI * 20.0 / 180.0);

    let dubins = DubinsManeuver3D::new(qi, qf, min_turn_radius, pitch_lims);
    let samples = dubins.compute_sampling(500);
}
```