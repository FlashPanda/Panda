#pragma once
#include "math/PandaMath.hpp"
#include "Filter.hpp"
#include "Spectrum.hpp"
#include <mutex>

namespace Panda
{
    struct FilmTilePixel
    {
        Spectrum contribSum = 0.f;
        float filterWeightSum = 0.f;
    };

    class FilmTile;

    class Film
    {
        public:
            Film(const Point2Di& resolution, const Bounds2Df& cropWindow,
                 std::unique_ptr<Filter> filter, float diagonal,
                 const std::string& filename, float scale,
                 float maxSampleLuminance = Infinity);
            Bounds2Di GetSampleBounds() const;
            Bounds2Df GetPhysicalExtent() const;
            std::unique_ptr<FilmTile> GetFilmTile(const Bounds2Di& sampleBounds);
            void MergeFilmTile(std::unique_ptr<FilmTile> tile);
            void SetImage(const Spectrum* img) const;
            void AddSplat(const Point2Df& p, Spectrum v);
            void WriteImage(float splatScale = 1);
            void Clear();

            const Point2Di m_FullResolution;
            const float m_Diagonal;
            std::unique_ptr<Filter> m_Filter;
            const std::string m_Filename;
            Bounds2Di m_CroppedPixelBounds;

        private:
            struct Pixel
            {
                Pixel() {xyz[0] = xyz[1] = xyz[2] = filterWeightSum = 0.f;}
                float xyz[3];
                float filterWeightSum;
                AtomicFloat splatXYZ[3];
                float pad;
            };
            std::unique_ptr<Pixel[]> pixels;
            static constexpr int32_t m_FilterTableWidth = 16;
            float m_FilterTable[m_FilterTableWidth * m_FilterTableWidth];
            std::mutex m_Mutex;
            const float m_Scale;
            const float m_MaxSampleLuminance;

            Pixel& GetPixel(const Point2Di& p)
            {
                assert(InsideExclusive(p, m_CroppedPixelBounds));
                int32_t width = m_CroppedPixelBounds.pMax[0] - m_CroppedPixelBounds.pMin[0];
                int32_t offset = (p[0] - m_CroppedPixelBounds.pMin[0]) +
                                 (p[1] - m_CroppedPixelBounds.pMin[1]) * width;
                return pixels[offset];
            }
    };

    class FilmTile
    {

    };
}