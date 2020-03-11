#pragma once
#include <utility>
#include <tiledarray.h>

template<typename Tile, typename Builder>
struct DirectTile {
    using eval_type = Tile;
    using numeric_type = typename Tile::numeric_type;

    TA::Range range;
    Builder builder;

    DirectTile() = default;
    DirectTile(const DirectTile& other) = default;
    DirectTile& operator=(const DirectTile& other) = default;
    DirectTile(TA::Range& range, Builder& builder) : range(std::move(range)), builder(std::move(builder)){}
    
    explicit operator eval_type() const {
        return builder(range);
    }
    
    template<typename Archive>
    void serialize(Archive &ar) {
        ar &range;
        ar &builder;
    }

}; // class DirectTile