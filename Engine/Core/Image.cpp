#include "Image.hpp"

namespace Panda
{
    std::ostream& operator<<(std::ostream& out, const Image& image)
    {
        out << "Image" << std::endl;
        out << "-----" << std::endl;
        out << "Width: " << image.Width << std::endl;
        out << "Height: " << image.Height << std::endl;
        out << "Bit Count: " << image.BitCount << std::endl;
        out << "Pitch: " << image.Pitch << std::endl;
        out << "Data Size: " << image.DataSize << std::endl;
#if DUMP_DETAILS
        int byteCount = image.BitCount >> 3;

        for (uint32_t i = 0; i < image.Height; ++i)
        {
            for (uint32_t j = 0; j < image.Width; ++j)
            {
                for (auto k = 0; k < byteCount; ++k)
                {
                    printf("%x ", reinterpret_cast<uint8_t*>(image.Data)[image.Pitch * i + j * byteCount + k]);   
                }
                std::cout << "\t";
            }
            std::cout << std::endl;
        }
#endif
        return out;
    }
}