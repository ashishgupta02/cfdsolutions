#define VERDICT_EXPORTS

#include "verdict.h"
#include <math.h>

/*!\brief Length of and edge.
 * Length is calculated by taking the distance between the end nodes.
 */
C_FUNC_DEF double v_edge_length(int /*num_nodes*/, double coordinates[][3]) {

    double x = coordinates[1][0] - coordinates[0][0];
    double y = coordinates[1][1] - coordinates[0][1];
    double z = coordinates[1][2] - coordinates[0][2];
    return (double) (sqrt(x * x + y * y + z * z));
}

/*!\brief Higher order function for calculating multiple metrics at once.
 * For an edge, there is only one metric, edge length.
 */
C_FUNC_DEF void v_edge_quality(int num_nodes, double coordinates[][3],
        unsigned int metrics_request_flag, struct EdgeMetricVals *metric_vals) {
    if (metrics_request_flag & V_EDGE_LENGTH)
        metric_vals->length = v_edge_length(num_nodes, coordinates);
}


