#include "Film.hpp"

namespace Panda
{
    Film::Film(const Point2Di& resolution, const Bounds2Df& cropWindow,
            std::unique_ptr<Filter> filter, float diagonal,
            const std::string& filename, float scale,
            float maxSampleLuminance)
        : m_FullResolution(resolution),  m_Diagonal(diagonal * 0.001f),
          m_Filter(std::move(filter)), m_Filename(filename),
          m_Scale(scale), m_MaxSampleLuminance(maxSampleLuminance)
    {}

    Bounds2Di Film::GetSampleBounds() const
    {
        return Bounds2Di();
    }

    Bounds2Df Film::GetPhysicalExtent() const
    {
		return Bounds2Df();
	}

    std::unique_ptr<FilmTile> Film::GetFilmTile(const Bounds2Di& sampleBounds)
    {
        return std::unique_ptr<FilmTile>(new FilmTile());
    }
    void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile)
    {}

    void Film::SetImage(const Spectrum* img) const
    {}

    void Film::AddSplat(const Point2Df& p, Spectrum v)
    {

    }

    void Film::WriteImage(float splatScale)
    {}

    void Film::Clear()
    {}
}