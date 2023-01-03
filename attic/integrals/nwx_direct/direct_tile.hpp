/*
 * Copyright 2023 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * Copyright 2022 NWChemEx-Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once
#include <memory>
#include <tiledarray.h>
#include <utility>

// A nearly general TiledArray lazy tile for use in direct methods.
template<typename Tile, typename Builder>
struct DirectTile {
    using eval_type    = Tile;
    using numeric_type = typename Tile::numeric_type;
    using value_type   = typename Tile::value_type;

    TA::Range range; // The range of the tile
    Builder builder; // The builder that produces the tile data on call

    DirectTile()                        = default;
    DirectTile(const DirectTile& other) = default;
    DirectTile& operator=(const DirectTile& other) = default;
    DirectTile(TA::Range& range, Builder builder) :
      range(std::move(range)), builder(std::move(builder)) {}

    /** @brief Convert to data tile type
     *
     *  @returns The filled data tile
     */
    explicit operator eval_type() { return builder(range); }

    /** @brief Serialize the tile
     *
     *  @param ar The archive
     */
    template<typename Archive>
    void serialize(Archive& ar) {
        ar& range;
        ar& builder;
    }

}; // class DirectTile

/** @brief Stream operator for Direct Tile
 *
 *  @param os Stream
 *  @param t The direct tile
 *  @returns Output Stream
 */
template<typename Tile, typename Builder>
std::ostream& operator<<(std::ostream& os, const DirectTile<Tile, Builder>& t) {
    os << t.range << "\n";
    return os;
}