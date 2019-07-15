
Mesh simplification implementation

- Restriction on input:
    - No non-manifold edge: will throw an exception
    - No non-manifold vertex: undefined behavior

- To-do-next:
    1. Use OOP better to improve program
        - neighborhood


Naming
    class:          Erasable
    typedef:        Indices
    member:         Faces.replaceSide(1)
    private member: _indices
    file name:      lowercase.{hpp,cpp}

Terms
    faces incident to some edge are its wings
    edges around a face are its sides
    vertices of an edge are its endpoints
    vertices of a face are its corners

=====
How was topology changed during simplification?
=====

After initial edge collapse operation on edge (d, k), two endpoints merge into
one, say k. Here it is possible that (a, d) and (a, k) two edges become coincide
which will introduce non-manifold edges in output.

This implementation eliminate non-manifold edges immediately after such
coincide edges appear. The method is *forking*, which separate the two edges
by attaching the forks of endpoints to one of the two coincide edges.

E.g., two edges ek, ed both have endpoints (k, a) but they have different wings.
First create new vertices k' and a' as forks. One of ed and ek will become
(k', a') instead of (k, a). Second, exchange one wing between the ed and ek
to consolidate the connectivity.
