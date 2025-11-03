# SafeMotionPlanner-FixedWingUAS
Supplementary material for “Development and Application of a Dynamic Obstacle Avoidance Algorithm for Small Fixed-Wing Aircraft with Safety Guarantees” (Control Engineering Practice, 2025). Includes appendix, code, and flight test data.

Abstract: This paper presents a motion planner for a small fixed-wing uncrewed aircraft system navigating in a dynamic obstacle environment. This planner generates dynamically feasible, evasive trajectories by concatenating motion primitives from a predefined library. A key contribution is the use of the dReal4 satisfiability modulo theories solver, which guarantees the existence of a safe evasive trajectory under a set of environmental assumptions. To reduce tracking errors due to the wind disturbance, an adaptive two-state extended Kalman filter is used to estimate steady wind and accordingly adjust the prespecified motion primitives. The controller used for implementation is a switched $\ell_2$-induced norm controller. The effectiveness of the proposed approach is validated through high-fidelity simulations and flight tests. The results from both simulation and flight tests are consistent with the safety guarantees obtained for the motion planner.

### Citation

If you use this repository, please cite as:

@article{Marquis2025,
  title={Development and Application of a Dynamic Obstacle Avoidance Algorithm for Small Fixed-Wing Aircraft with Safety Guarantees},
  author={Dennis J. Marquis and Mazen Farhood},
  journal={Control Engineering Practice},
  year={2025},
  note={Code and supplementary materials available at \url{https://github.com/dennisjmarquis/SafeMotionPlanner-FixedWingUAS}}
}


### License

This repository is licensed under the MIT License. See `LICENSE` for details.
