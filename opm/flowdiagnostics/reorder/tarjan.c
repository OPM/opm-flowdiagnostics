/*
Copyright (C) 2012 (c) Jostein R. Natvig <jostein natvig at gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <opm/flowdiagnostics/reorder/tarjan.h>

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

enum VertexMark { DONE = -2, REMAINING = -1 };

struct TarjanWorkSpace
{
    int nvert;

    int *status;
    int *link;
    int *time;

    int *stack;
    int *cstack;
};

static void
initialise_stacks(int *comp, int *vert, struct TarjanWorkSpace *ws)
{
    assert (comp != vert);

    ws->stack  = comp + (ws->nvert + 0);
    ws->cstack = vert + (ws->nvert - 1);
}

static void
assign_int_vector(size_t n, const int val, int *v)
{
    size_t i;

    for (i = 0; i < n; i++) { v[i] = val; }
}

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* Stack grows to lower addresses */
#define peek(stack) (*((stack) + 1))
#define push(stack) *(stack)--

/* ======================================================================
 * Public interface below separator
 * ====================================================================== */

struct TarjanWorkSpace *
create_tarjan_workspace(const int nvert)
{
    struct TarjanWorkSpace *ws, ws0 = { 0 };

    ws = malloc(1 * sizeof *ws);

    if (ws != NULL) {
        *ws = ws0;

        ws->status = malloc(3 * nvert * sizeof *ws->status);

        if (ws->status == NULL) {
            destroy_tarjan_workspace(ws);

            ws = NULL;
        }
        else {
            ws->nvert = nvert;

            ws->link = ws->status + ws->nvert;
            ws->time = ws->link   + ws->nvert;
        }
    }

    return ws;
}

void
destroy_tarjan_workspace(struct TarjanWorkSpace *ws)
{
    if (ws != NULL) {
        free(ws->status);
    }

    free(ws);
}

/*
  Compute the strong components of directed graph G(edges, vertices),
  return components in reverse topological sorted sequence.
  Complexity O(|vertices|+|edges|). See "http://en.wikipedia.org/wiki/
  Tarjan's_strongly_connected_components_algorithm".

  nv    - number of vertices

  ia,ja - adjacency matrix for directed graph in compressed sparse row
          format: vertex i has directed edges to vertices ja[ia[i]],
          ..., ja[ia[i+1]-1].

  vert  - permutation of vertices into topologically sorted sequence of
          strong components (i.e., loops).

  comp  - pointers to start of each strongly connected component in
          vert, the i'th component has vertices vert[comp[i]], ...,
          vert[comp[i+1]-1].

  ncomp - number of strong components.

  work  - block of memory of size 3*nv*sizeof(int).
 */

/*--------------------------------------------------------------------*/
void
tarjan(int                     nv,
       const int              *ia,
       const int              *ja,
       int                    *vert,
       int                    *comp,
       int                    *ncomp,
       struct TarjanWorkSpace *work)
/*--------------------------------------------------------------------*/
{
    int  c, v, seed, child, t, pos;

    /*
     * Note: 'status' serves dual purpose during processing.
     *
     *   status[c] = DONE      -> Cell 'c' fully processed.
     *               REMAINING -> Cell 'c' Undiscovered.
     *               0         -> Cell 'c' not fully classified--there are
     *                            no remaining descendants for this cell but
     *                            we do not yet know if the cell is an SCC
     *                            all by itself or if it is part of a larger
     *                            component.
     *
     *   status[c] > 0 -> Cell 'c' has 'status[c]' remaining
     *                    descendants (i.e., successors/children)
     */

    int *status = NULL;
    int *link   = NULL;
    int *time   = NULL;
    int *stack  = NULL;
    int *cstack = NULL;

    if (nv != work->nvert) {
        return;
    }

    /* Hint: end of VERT and COMP are used as stacks. */

    initialise_stacks(comp, vert, work);

    stack  = work->stack;
    cstack = work->cstack;

    status = work->status;
    link   = work->link;
    time   = work->time;

    /* Init status all vertices */
    assign_int_vector(nv, REMAINING, status);

    *ncomp  = 0;
    *comp++ = pos = 0;

    seed = 0;
    while (seed < nv)
    {
        if (status[seed] == DONE) {
            ++seed;
            continue;
        }

        push(stack) = seed;

        t = 0;

        while (stack != work->stack)
        {
            c = peek(stack);

            assert(status[c] != DONE);
            assert(status[c] >= -2);

            if (status[c] == REMAINING) {
                /* number of descendants of c */
                status[c] = ia[c + 1] - ia[c];
                time[c]   = link[c] = t++;

                push(cstack) = c;
            }

            /* if all descendants are processed */
            if (status[c] == 0)
            {
                /* if c is root of strong component */
                if (link[c] == time[c])
                {
                    do {
                        assert (cstack != work->cstack);

                        /* pop strong component stack */
                        v         = *++cstack;
                        status[v] = DONE;

                        /* store vertex in VERT */
                        vert[pos++]  = v;
                    } while (v != c);

                    /* store end point of component */
                    *comp++ = pos;
                    *ncomp += 1;
                }

                /* pop c */
                ++stack;

                if (stack != work->stack) {
                    v = peek(stack);

                    link[v] = MIN(link[v], link[c]);
                }
            }

            /* if there are more descendants to consider */
            else {
                assert(status[c] > 0);

                child = ja[ia[c] + (status[c] - 1)];

                /* decrement descendant count of c*/
                --status[c];

                if (status[child] == REMAINING) {
                    /* push child */
                    push(stack) = child;
                }
                else if (status[child] >= 0) {
                    link[c] = MIN(link[c], time[child]);
                }
                else {
                    assert(status[child] == DONE);
                }
            }
        }

        assert (cstack == work->cstack);
    }
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
