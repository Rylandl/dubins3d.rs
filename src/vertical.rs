use core::cmp::Ordering;

use crate::{PI, mod2pi};
use crate::dubins2d::{
    DubinsManeuver2D,
    DubinsStruct,
    ManeuverCase,
    SegmentType
};

pub(crate) fn get_vertical(qi: (f64, f64, f64), qf: (f64, f64, f64), rhomin: f64, pitchmax: (f64, f64)) -> DubinsManeuver2D {
    let mut maneuver = DubinsManeuver2D{
        qi: qi, 
        qf: qf,
        rhomin: rhomin,
        maneuver: DubinsStruct{
            t: 0.0,
            p: 0.0,
            q: 0.0, 
            length: core::f64::INFINITY,
            case: ManeuverCase{a: SegmentType::NONE, b: SegmentType::NONE, c: SegmentType::NONE}
        }
    };

    let path_lsl = _lsl(&maneuver);
    let path_rsr = _rsr(&maneuver);
    let path_lsr = _lsr(&maneuver, pitchmax);
    let path_rsl = _rsl(&maneuver, pitchmax);
    let mut _paths = vec![path_lsr, path_lsl, path_rsr, path_rsl];
    
    // _paths.sort_by(|a,b| a.length.partial_cmp(&b.length).unwrap());
    _paths.sort_by(|p1, p2| {
        let a = p1.length;
        let b = p2.length;
        match (a.is_nan(), b.is_nan()) {
            (true, true) => Ordering::Equal,
            (true, false) => Ordering::Greater,
            (false, true) => Ordering::Less,
            (false, false) => a.partial_cmp(&b).unwrap(),
        }
    });

    for path in _paths {
        if path.t.abs() < PI && path.q.abs() < PI {
            let center_angle: f64;
            if path.case.a == SegmentType::LEFT {
                center_angle = maneuver.qi.2 + path.t;
            }
            else {
                center_angle = maneuver.qi.2 - path.t;
            }
            if center_angle < pitchmax.0 || center_angle > pitchmax.1 {
                continue;
            }
            maneuver.maneuver = path;
            break;
        }
    }

    return maneuver;
}

fn _lsl(maneuver: &DubinsManeuver2D) -> DubinsStruct {
    let theta1 = maneuver.qi.2;
    let theta2 = maneuver.qf.2;

    let mut t = core::f64::INFINITY;
    let mut p = core::f64::INFINITY;
    let mut q = core::f64::INFINITY;
    if theta1 <= theta2 {
        let p1 = (maneuver.qi.0, maneuver.qi.1);
        let p2 = (maneuver.qf.0, maneuver.qf.1);

        let radius = maneuver.rhomin;

        let c1 = radius * theta1.cos();
        let s1 = radius * theta1.sin();
        let c2 = radius * theta2.cos();
        let s2 = radius * theta2.sin();

        let o1 = (p1.0 - s1, p1.1 + c1);
        let o2 = (p2.0 - s2, p2.1 + c2);

        let diff = (o2.0 - o1.0, o2.1 - o1.1);
        let center_distance = (diff.0 * diff.0 + diff.1 * diff.1).sqrt();
        let center_angle = (diff.1).atan2(diff.0);

        t = mod2pi(-theta1 + center_angle);
        p = center_distance / radius;
        q = mod2pi(theta2 - center_angle);

        if t > PI {
            t = 0.0;
            q = theta2 - theta1;
            let turn_end_y = o2.1 - radius * theta1.cos();
            let diff_y = turn_end_y - p1.1;
            if theta1.abs() > 1e-5 && ((diff_y < 0.0) == (theta1 < 0.0)) {
                p = diff_y / theta1.sin() / radius;
            }
            else {
                t = core::f64::INFINITY;
                p = core::f64::INFINITY;
                q = core::f64::INFINITY;
            }
        }
        if q > PI {
            t = theta2 - theta1;
            q = 0.0;
            let turn_end_y = o1.1 - radius * theta2.cos();
            let diff_y = p2.1 - turn_end_y;
            if theta2.abs() > 1e-5 && ((diff_y < 0.0) == (theta2 < 0.0)) {
                p = diff_y / theta2.sin() / radius;
            }
            else {
                t = core::f64::INFINITY;
                p = core::f64::INFINITY;
                q = core::f64::INFINITY;
            }
        }
    }

    let length = (t+p+q) * maneuver.rhomin;
    
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {
            a: SegmentType::LEFT,
            b: SegmentType::STRAIGHT,
            c: SegmentType::LEFT
        }
    };
    return ds;
}

fn _rsr(maneuver: &DubinsManeuver2D) -> DubinsStruct {
    let theta1 = maneuver.qi.2;
    let theta2 = maneuver.qf.2;

    let mut t = core::f64::INFINITY;
    let mut p = core::f64::INFINITY;
    let mut q = core::f64::INFINITY;
    if theta2 <= theta1 {
        let p1 = (maneuver.qi.0, maneuver.qi.1);
        let p2 = (maneuver.qf.0, maneuver.qf.1);

        let radius = maneuver.rhomin;

        let c1 = radius * theta1.cos();
        let s1 = radius * theta1.sin();
        let c2 = radius * theta2.cos();
        let s2 = radius * theta2.sin();

        let o1 = (p1.0 + s1, p1.1 - c1);
        let o2 = (p2.0 + s2, p2.1 - c2);

        let diff = (o2.0 - o1.0, o2.1 - o1.1);
        let center_distance = (diff.0 * diff.0 + diff.1 * diff.1).sqrt();
        let center_angle = (diff.1).atan2(diff.0);

        t = mod2pi(theta1 - center_angle);
        p = center_distance / radius;
        q = mod2pi(-theta2 + center_angle);

        if t > PI {
            t = 0.0;
            q = -theta2 + theta1;
            let turn_end_y = o2.1 + radius * theta1.cos();
            let diff_y = turn_end_y - p1.1;
            if theta1.abs() > 1e-5 && (diff_y < 0.0) == (theta1 < 0.0) {
                p = diff_y / theta1.sin() / radius;
            }
            else {
                t = core::f64::INFINITY;
                p = core::f64::INFINITY;
                q = core::f64::INFINITY;
            }
        }
        if q > PI {
            t = -theta2 + theta1;
            q = 0.0;
            let turn_end_y = o1.1 + radius * theta2.cos();
            let diff_y = p2.1 - turn_end_y;
            if theta2.abs() > 1e-5 && (diff_y < 0.0) == (theta2 < 0.0) {
                p = diff_y / theta2.sin() / radius;
            }
            else {
                t = core::f64::INFINITY;
                p = core::f64::INFINITY;
                q = core::f64::INFINITY;
            }
        }
    }

    let length = (t+p+q) * maneuver.rhomin;
    
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {
            a: SegmentType::RIGHT,
            b: SegmentType::STRAIGHT,
            c: SegmentType::RIGHT
        }
    };
    return ds;
}


fn _lsr(maneuver: &DubinsManeuver2D, pitchmax: (f64, f64)) -> DubinsStruct {
    let theta1 = maneuver.qi.2;
    let theta2 = maneuver.qf.2;

    let t: f64;
    let p: f64;
    let q: f64;

    let p1 = (maneuver.qi.0, maneuver.qi.1);
    let p2 = (maneuver.qf.0, maneuver.qf.1);

    let radius = maneuver.rhomin;

    let c1 = radius * theta1.cos();
    let s1 = radius * theta1.sin();
    let c2 = radius * theta2.cos();
    let s2 = radius * theta2.sin();

    let o1 = (p1.0 - s1, p1.1 + c1);
    let o2 = (p2.0 + s2, p2.1 - c2);

    let mut diff = (o2.0 - o1.0, o2.1 - o1.1);
    let center_distance = (diff.0 * diff.0 + diff.1 * diff.1).sqrt();
    
    let mut alpha = (2.0 * radius / center_distance).asin();
    if center_distance < 2.0 * radius {
        diff.0 = (4.0 * radius * radius - diff.1 *diff.1).sqrt();
        alpha = PI/2.0;
    }
    
    let mut center_angle = (diff.1).atan2(diff.0) + alpha;

    if center_angle < pitchmax.1 {
        t = mod2pi(-theta1 + center_angle);
        p = (center_distance * center_distance - 4.0 * radius * radius).max(0.0).sqrt() / radius;
        q = mod2pi(-theta2 + center_angle);
    }
    else {
        center_angle = pitchmax.1;
        t = mod2pi(-theta1 + center_angle);
        q = mod2pi(-theta2 + center_angle);

        let c = radius * center_angle.cos();
        let s = radius * center_angle.sin();
        let w1 = (o1.0 + s, o1.1 - c);
        let w2 = (o2.0 - s, o2.1 + c);

        p = (w2.1 - w1.1) / center_angle.sin() / radius;
    }

    let length = (t+p+q) * maneuver.rhomin;
    
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {
            a: SegmentType::LEFT,
            b: SegmentType::STRAIGHT,
            c: SegmentType::RIGHT
        }
    };
    return ds;
}

fn _rsl(maneuver: &DubinsManeuver2D, pitchmax: (f64, f64)) -> DubinsStruct {
    let theta1 = maneuver.qi.2;
    let theta2 = maneuver.qf.2;

    let t: f64;
    let p: f64;
    let q: f64;

    let p1 = (maneuver.qi.0, maneuver.qi.1);
    let p2 = (maneuver.qf.0, maneuver.qf.1);

    let radius = maneuver.rhomin;

    let c1 = radius * theta1.cos();
    let s1 = radius * theta1.sin();
    let c2 = radius * theta2.cos();
    let s2 = radius * theta2.sin();

    let o1 = (p1.0 + s1, p1.1 - c1);
    let o2 = (p2.0 - s2, p2.1 + c2);

    let mut diff = (o2.0 - o1.0, o2.1 - o1.1);
    let center_distance = (diff.0 * diff.0 + diff.1 * diff.1).sqrt();
    
    let mut alpha = (2.0 * radius / center_distance).asin();
    if center_distance < 2.0 * radius {
        diff.0 = (4.0 * radius * radius - diff.1 * diff.1).sqrt();
        alpha = PI/2.0;
    }
    
    let mut center_angle = (diff.1).atan2(diff.0) - alpha;

    if center_angle > pitchmax.0 {
        t = mod2pi(theta1 - center_angle);
        p = (center_distance * center_distance - 4.0 * radius * radius).max(0.0).sqrt() / radius;
        q = mod2pi(theta2 - center_angle);
    }
    else {
        center_angle = pitchmax.0;
        t = mod2pi(theta1 - center_angle);
        q = mod2pi(theta2 - center_angle);

        let c = radius * center_angle.cos();
        let s = radius * center_angle.sin();
        let w1 = (o1.0 - s, o1.1 + c);
        let w2 = (o2.0 + s, o2.1 - c);

        p = (w2.1 - w1.1) / center_angle.sin() / radius;
    }

    let length = (t+p+q) * maneuver.rhomin;
    
    let ds = DubinsStruct {
        t: t,
        p: p,
        q: q,
        length: length,
        case: ManeuverCase {
            a: SegmentType::RIGHT,
            b: SegmentType::STRAIGHT,
            c: SegmentType::LEFT
        }
    };
    return ds;
}
