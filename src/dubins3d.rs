#![allow(dead_code)]

// mod dubins2d;
use crate::dubins2d::{
    DubinsManeuver2D, 
    SegmentType, 
    ManeuverCase,
    get_coordinates_at
};
use crate::vertical;

// use dubins2d::DubinsManeuver2D

type Point = (f64, f64, f64, f64, f64);

pub struct DubinsManeuver3D {
    qi: Point,
    qf: Point,
    rhomin: f64,
    pitchlims: (f64, f64),
    path: Vec<DubinsManeuver2D>,
    length: f64
}

impl DubinsManeuver3D {
    pub fn new(qi: Point, qf: Point, rhomin: f64, pitchlims: (f64, f64)) -> DubinsManeuver3D {
        let mut maneuver = DubinsManeuver3D {
            qi: qi,
            qf: qf,
            rhomin: rhomin,
            pitchlims: pitchlims,
            path: Vec::new(),
            length: -1.0
        };
    
        let a = 1.0;
        let mut b = 1.0;
    
        let fa = try_to_construct(&maneuver, maneuver.rhomin * a);
        let mut fb = try_to_construct(&maneuver, maneuver.rhomin * b);
    
        while fb.len() < 2 {
            b *= 2.0;
            fb = try_to_construct(&maneuver, maneuver.rhomin * b);
        }
    
        if fa.len() > 0 {
            maneuver.path = fa;
        }
        else {
            if fb.len() < 2 {
                // return Err(
                //     Error::new(ErrorKind::)
                // );
                return maneuver;
            }
        }
    
        let mut step: f64 = 0.1;
        while step.abs() > 1e-10 {
            let c = (b + step).max(1.0);
            let fc = try_to_construct(&maneuver, maneuver.rhomin * c);
            if fc.len() > 0 {
                if fc[1].maneuver.length < fb[1].maneuver.length {
                    b = c;
                    fb = fc;
                    step *= 2.0;
                    continue;
                }
            }
            step *= -0.1;
        }
        maneuver.path.clear();
        maneuver.path.extend(fb);
        maneuver.length = maneuver.path[1].maneuver.length;
    
        return maneuver;
    }
    
    pub fn get_lower_bound(qi: Point, qf: Point, rhomin: f64, pitchlims: (f64, f64)) -> Self {
        let mut maneuver = DubinsManeuver3D {
            qi: qi,
            qf: qf,
            rhomin: rhomin,
            pitchlims: pitchlims,
            path: Vec::new(),
            length: -1.0
        };
    
        let spiral_radius = rhomin * ((-pitchlims.0).max(pitchlims.1)).cos().powi(2);
    
        let qi2d = (maneuver.qi.0, maneuver.qi.1, maneuver.qi.3);
        let qf2d = (maneuver.qf.0, maneuver.qf.1, maneuver.qf.3);
        let dlat = DubinsManeuver2D::new(qi2d, qf2d, spiral_radius, core::f64::NEG_INFINITY, false);
    
        let qi3d = (0.0, maneuver.qi.2, maneuver.qi.4);
        let qf3d = (dlat.maneuver.length, maneuver.qf.2, maneuver.qf.4);
        let dlon = vertical::get_vertical(qi3d, qf3d, maneuver.rhomin, maneuver.pitchlims);
    
        if dlon.maneuver.case.a == SegmentType::NONE {
            maneuver.length = 0.0;
            return maneuver;
        }
    
        maneuver.length = dlon.maneuver.length;
        maneuver.path.extend([dlat, dlon]);
        return maneuver;
    }
    
    pub fn get_upper_bound(qi: Point, qf: Point, rhomin: f64, pitchlims: (f64, f64)) -> DubinsManeuver3D {
        let mut maneuver = DubinsManeuver3D {
            qi: qi,
            qf: qf,
            rhomin: rhomin,
            pitchlims: pitchlims,
            path: Vec::new(),
            length: -1.0
        };
    
        let safe_radius = (2.0 as f64).sqrt() * maneuver.rhomin;
        let diff = (qf.0 - qi.0, qf.1 - qi.1);
        let dist = (diff.0 * diff.0 + diff.1 * diff.1).sqrt();
        if dist < 4.0 * safe_radius {
            maneuver.length = core::f64::INFINITY;
            return maneuver;
        }
    
        let qi2d = (maneuver.qi.0, maneuver.qi.1, maneuver.qi.3);
        let qf2d = (maneuver.qf.0, maneuver.qf.1, maneuver.qf.3);
        let dlat = DubinsManeuver2D::new(qi2d, qf2d, safe_radius, core::f64::NEG_INFINITY, false);
    
        let qi3d = (0.0, maneuver.qi.2, maneuver.qi.4);
        let qf3d = (dlat.maneuver.length, maneuver.qf.2, maneuver.qf.4);
        let dlon = vertical::get_vertical(qi3d, qf3d, safe_radius, maneuver.pitchlims);
    
        if dlon.maneuver.case.a == SegmentType::NONE {
            maneuver.length = core::f64::INFINITY;
            return maneuver;
        }
    
        maneuver.length = dlon.maneuver.length;
        maneuver.path.extend([dlat, dlon]);
        return maneuver;
    }

    pub fn compute_sampling(self: &Self, number_of_samples: i32) -> Vec<Point> {
        let dlat = &self.path[0];
        let dlon = &self.path[1];
    
        let mut points: Vec<Point> = Vec::new();
    
        for sample in 0..number_of_samples {
            let prog: f64 =  dlon.maneuver.length * (sample as f64) / (number_of_samples as f64);
            // println!("{}", prog);
            let q_sz = get_coordinates_at(&dlon, prog);
            let q_xy = get_coordinates_at(&dlat, q_sz.0);
            points.push((q_xy.0, q_xy.1, q_sz.1, q_xy.2, q_sz.2));
        }
    
        return points;
    }
}

fn try_to_construct(maneuver: &DubinsManeuver3D, horizontal_radius: f64) -> Vec<DubinsManeuver2D> {
    let qi2d = (maneuver.qi.0, maneuver.qi.1, maneuver.qi.3);
    let qf2d = (maneuver.qf.0, maneuver.qf.1, maneuver.qf.3);

    let dlat = DubinsManeuver2D::new(qi2d, qf2d, horizontal_radius, core::f64::NEG_INFINITY, false);

    let qi3d = (0.0, maneuver.qi.2, maneuver.qi.4);
    let qf3d = (dlat.maneuver.length, maneuver.qf.2, maneuver.qf.4);

    let vertical_curvature = (1.0 / maneuver.rhomin / maneuver.rhomin - 1.0 / horizontal_radius / horizontal_radius).sqrt();
    if vertical_curvature < 1e-5 {
        return vec![];
    }

    let vertical_radius = 1.0 / vertical_curvature;
    
    let dlon = DubinsManeuver2D::new(qi3d, qf3d, vertical_radius, core::f64::NEG_INFINITY, false);

    if dlon.maneuver.case == (ManeuverCase{a: SegmentType::RIGHT, b: SegmentType::LEFT, c: SegmentType::RIGHT}) ||
        dlon.maneuver.case == (ManeuverCase{a: SegmentType::RIGHT, b: SegmentType::LEFT, c: SegmentType::RIGHT}) {
        return vec![];
    }

    if dlon.maneuver.case.a == SegmentType::RIGHT {
        if maneuver.qi.4 - dlon.maneuver.t < maneuver.pitchlims.0 {
            return vec![];
        }
    }
    else {
        if maneuver.qi.4 + dlon.maneuver.t > maneuver.pitchlims.1 {
            return vec![];
        }
    }
    return vec![dlat, dlon];
}
