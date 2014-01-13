#include "pointTopo.h"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointTopo::pointTopo
(
    const std::set<std::set<Foam::label> > &triangles,
    blockMeshTopology *topo
)
    :
      triangles_(triangles),
      topo_(topo)
{
}

Foam::pointTopo::~pointTopo()
{

}

// * * * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
