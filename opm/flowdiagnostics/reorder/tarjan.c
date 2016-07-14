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
};

struct TarjanSCCResult
{
    size_t  ncomp;
    size_t *start;
    int    *vert;

    size_t *vstack, *vstart;
    int    *cstack, *cstart;
};

static void
initialise_stacks(const size_t nvert, struct TarjanSCCResult *scc)
{
    /* Vertex and component stacks grow from high end downwards */

    scc->vstack = scc->vstart = scc->start + (nvert + 0);
    scc->cstack = scc->cstart = scc->vert  + (nvert - 1);
}

static void
assign_int_vector(size_t n, const int val, int *v)
{
    size_t i;

    for (i = 0; i < n; i++) { v[i] = val; }
}

static struct TarjanSCCResult *
allocate_scc_result(const size_t nvert)
{
    struct TarjanSCCResult *scc, scc0 = { 0 };

    scc = malloc(1 * sizeof *scc);

    if (scc != NULL) {
        *scc = scc0;

        scc->start = malloc((nvert + 1) * sizeof *scc->start);
        scc->vert  = malloc(nvert       * sizeof *scc->vert);

        if ((scc->start == NULL) || (scc->vert == NULL)) {
            destroy_tarjan_sccresult(scc);

            scc = NULL;
        }
    }

    return scc;
}

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* Stack grows to lower addresses */
#define peek(stack) (*((stack) + 1))
#define push(stack) *(stack)--

static void
tarjan_global(const size_t            nv,
              const int              *ia,
              const int              *ja,
              struct TarjanWorkSpace *work,
              struct TarjanSCCResult *scc)
{
    int    t, vertex;
    size_t c, v, seed, child, pos;

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

    size_t *start  = NULL;

    assert ((work != NULL) && "Work array must be non-NULL");
    assert ((scc  != NULL) && "Result array must be non-NULL");

    initialise_stacks(work->nvert, scc);

    start  = scc->start;

    status = work->status;
    link   = work->link;
    time   = work->time;

    /* Init status all vertices */
    assign_int_vector(nv, REMAINING, status);

    scc->ncomp = 0;
    *start++   = pos = 0;

    seed = 0;
    while (seed < nv)
    {
        if (status[seed] == DONE) {
            ++seed;
            continue;
        }

        push(scc->vstack) = seed;

        t = 0;

        while (scc->vstack != scc->vstart)
        {
            c = peek(scc->vstack);

            assert(status[c] != DONE);
            assert(status[c] >= -2);

            if (status[c] == REMAINING) {
                /* number of descendants of c */
                status[c] = ia[c + 1] - ia[c];
                time[c]   = link[c] = t++;

                push(scc->cstack) = (int) c;
            }

            /* if all descendants are processed */
            if (status[c] == 0)
            {
                /* if c is root of strong component */
                if (link[c] == time[c])
                {
                    do {
                        assert (scc->cstack != scc->cstart);

                        /* pop strong component stack */
                        v = vertex = *++scc->cstack;
                        status[v]  = DONE;

                        /* store vertex in VERT */
                        scc->vert[pos++] = vertex;
                    } while (v != c);

                    /* store end point of component */
                    *start++    = pos;
                    scc->ncomp += 1;
                }

                /* pop c */
                ++scc->vstack;

                if (scc->vstack != scc->vstart) {
                    v = peek(scc->vstack);

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
                    push(scc->vstack) = child;
                }
                else if (status[child] >= 0) {
                    link[c] = MIN(link[c], time[child]);
                }
                else {
                    assert(status[child] == DONE);
                }
            }
        }

        assert (scc->cstack == scc->cstart);
    }
}

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

void
destroy_tarjan_sccresult(struct TarjanSCCResult *scc)
{
    if (scc != NULL) {
        /* Reverse order of acquisition. */
        free(scc->vert);
        free(scc->start);
    }

    free(scc);
}

size_t
tarjan_get_numcomponents(const struct TarjanSCCResult *scc)
{
    return scc->ncomp;
}

struct TarjanComponent
tarjan_get_strongcomponent(const struct TarjanSCCResult *scc,
                           const size_t                  compID)
{
    size_t start;

    assert ((compID < tarjan_get_numcomponents(scc)) &&
            "Component ID out of bounds");

    struct TarjanComponent c = { 0 };

    start = scc->start[compID];

    c.size   = scc->start[compID + 1] - start;
    c.vertex = &scc->vert[start];

    return c;
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
 */

/*--------------------------------------------------------------------*/
struct TarjanSCCResult *
tarjan(const int  nv,
       const int *ia,
       const int *ja)
/*--------------------------------------------------------------------*/
{
    struct TarjanWorkSpace *work = NULL;
    struct TarjanSCCResult *scc  = NULL;

    work = create_tarjan_workspace(nv);
    scc  = allocate_scc_result(nv);

    if ((work == NULL) || (scc == NULL)) {
        destroy_tarjan_sccresult(scc);
        destroy_tarjan_workspace(work);

        work = NULL;
        scc  = NULL;
    }

    if (scc != NULL) {
        tarjan_global(nv, ia, ja, work, scc);
    }

    destroy_tarjan_workspace(work);

    return scc;
}

/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
