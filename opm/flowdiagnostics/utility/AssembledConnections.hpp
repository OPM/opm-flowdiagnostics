/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_ASSEMBLEDCONNECTIONS_HEADER_INCLUDED
#define OPM_ASSEMBLEDCONNECTIONS_HEADER_INCLUDED

#include <vector>

namespace Opm { namespace Utility {

    class AssembledConnections
    {
    public:
        void addConnection(const int i, const int j);
        void addConnection(const int i, const int j, const double v);

        void compress();

        using Neighbours = std::vector<int>;
        using Offset     = Neighbours::size_type;
        using Start      = std::vector<Offset>;
        using ConnWeight = std::vector<double>;

        const Start& startPointers() const;

        const Neighbours& neighbourhood() const;

        const ConnWeight& connectionWeight() const;

    private:
        class Connections
        {
        public:
            using EntityVector = std::vector<int>;
            using WeightVector = std::vector<double>;

            void add(const int i, const int j);
            void add(const int i, const int j, const double v);

            void clear();

            bool empty() const;

            bool isValid() const;

            bool isWeighted() const;

            int maxRow() const;
            int maxCol() const;

            EntityVector::size_type nnz() const;

            const EntityVector& i() const;
            const EntityVector& j() const;
            const WeightVector& v() const;

        private:
            EntityVector i_;
            EntityVector j_;
            WeightVector v_;

            int max_i_{ -1 };
            int max_j_{ -1 };
        };

        class CSR
        {
        public:
            void create(const Connections& conns);

            const Start&      ia() const;
            const Neighbours& ja() const;
            const ConnWeight& sa() const;

        private:
            Start      ia_;
            Neighbours ja_;
            ConnWeight sa_;

            Start      elmIdx_;

            int        numRows_{ 0 };
            int        numCols_{ 0 };

            // ---------------------------------------------------------
            // Implementation of create()
            // ---------------------------------------------------------

            void assemble(const Connections& conns);

            void sort();

            void condenseDuplicates();

            void accumulateConnWeights(const std::vector<double>& v);

            // ---------------------------------------------------------
            // Implementation of assemble()
            // ---------------------------------------------------------

            void accumulateRowEntries(const int               numRows,
                                      const std::vector<int>& rowIdx);

            void createGraph(const std::vector<int>& rowIdx,
                             const std::vector<int>& colIdx);

            // ---------------------------------------------------------
            // General utilities
            // ---------------------------------------------------------

            void transpose();

            std::vector<int> expandStartPointers() const;

            void unique(Neighbours::const_iterator begin,
                        Neighbours::const_iterator end);

            void remapElementIndex(Start&& elmIdx);
        };

        Connections conns_;
        CSR         csr_;
    };

}} // namespace Opm::Utility

#endif // OPM_ASSEMBLEDCONNECTIONS_HEADER_INCLUDED