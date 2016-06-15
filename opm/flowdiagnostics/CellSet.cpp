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

#if HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <opm/flowdiagnostics/CellSet.hpp>

#include <utility>

Opm::CellSetID::CellSetID()
{}

Opm::CellSetID::CellSetID(Repr id)
    : id_(std::move(id))
{
}

std::string
Opm::CellSetID::to_string() const
{
    return id_;
}

// =====================================================================

void
Opm::CellSet::identify(CellSetID id)
{
    id_ = std::move(id);
}

const Opm::CellSetID&
Opm::CellSet::id() const
{
    return id_;
}

void
Opm::CellSet::insert(const int i)
{
    iset_.insert(i);
}

Opm::CellSet::const_iterator
Opm::CellSet::begin() const
{
    return iset_.begin();
}

Opm::CellSet::const_iterator
Opm::CellSet::end() const
{
    return iset_.end();
}
