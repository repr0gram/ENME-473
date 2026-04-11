 (or as close as possible). That's geometrically intuitive — if you can put Pin A anywhere on link 8, the           
    box-minimizing choice is to put it where it doesn't protrude. The 50% reduction is real in the sense that the math  
    is exact, but it's not a meaningful design: it corresponds to a link-8 triangle that's essentially a straight line, 
    which probably isn't physically realizable (you can't pin three points on a rigid body if two of them are the same  
    point, and a near-degenerate triangle has terrible stiffness).                                                      
   
    What you should probably do:                                                                                        
                                                            
    - Add a minimum triangle "fatness" constraint, e.g. require |pinA_yL| ≥ some_min (say 50 mm or 100 mm), so Pin A has
     meaningful perpendicular offset from the P1–P2 axis. This reflects the actual design reality that link 8 needs to
    be a real triangle with a real Pin A sticking off it.
    - Or constrain the minimum interior angle of the P1–P2–PinA triangle away from 0° and 180°.
    - Re-run and see what the best design looks like with, say, |pinA_yL| ≥ 100 mm. My guess is you'll get something in
    the 25–35% area-reduction range — worse than 50% but actually usable, and probably similar in character to the
    forward search's result (which was constrained implicitly because its grid didn't explore the degenerate edge).

    TL;DR: the 50% number is real but the geometry it corresponds to is degenerate. The forward search's 24% answer is
    almost certainly closer to what a realistic constrained optimum looks like. Adding a fatness constraint on the
    link-8 triangle will give you a meaningful answer that sits between those two.
