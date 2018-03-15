
 OPCODE distribution 1.0
 -----------------------

 WHAT ?

    OPCODE means OPtimized COllision DEtection.
    So this is a collision detection package similar to RAPID. Here's a
    quick list of features:

    - C++ interface, developed for Windows systems using VC++ 6.0
    - Works on arbitrary meshes (convex or non-convex), even polygon soups
    - Current implementation uses AABB-trees
    - Supports CD queries for 2 bodies only. The sweep-and-prune and N-body code is part
      of a bigger unreleased package (Z-COLLIDE)
    - Introduces Primitive-BV overlap tests during recursive collision queries (whereas
      standard libraries only rely on Primitive-Primitive and BV-BV tests)
    - Introduces no-leaf trees, i.e. collision trees whose leaf nodes have been removed
    - Supports collision queries on quantized trees (decompressed on-the-fly)
    - Supports "first contact" or "all contacts" modes (à la RAPID)
    - Uses temporal coherence for "first contact" mode (~10 to 20 times faster, useful
      in rigid body simulation during bisection)
    - Memory footprint is 7.2 times smaller than RAPID's one, which is ideal for console
      games with limited ram (actually, if you use the unmodified RAPID code using double
      precision, it's more like 13 times smaller...)
    - And yet it often runs faster than RAPID (according to RDTSC, sometimes more than 5
      times faster when objects are deeply overlapping)
    - Performance is usually close to RAPID's one in close-proximity situations
    - Does not work on deformable meshes. Yet.

 WHY ?

    - Because RAPID uses too many bytes.
    - Because the idea was nice...

 WHEN ?

    It's been coded in march 2001 following a thread on the GD-Algorithms list.

      GDAlgorithms-list mailing list
      GDAlgorithms-list@lists.sourceforge.net
      http://lists.sourceforge.net/lists/listinfo/gdalgorithms-list

 WHO ?

    Pierre Terdiman
    May 03, 2001

    p.terdiman@wanadoo.fr
    p.terdiman@codercorner.com
 
    http://www.codercorner.com
    http://www.codercorner.com/Opcode.htm
