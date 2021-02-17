#pragma once
#include "BaseSceneObject.hpp"
#include "Core/Math/PandaMath.hpp"
#include "Core/Parser/JPEG.hpp"
#include "Core/Parser/PNG.hpp"
#include "Core/Parser/BMP.hpp"
#include "Core/Parser/TGA.hpp"
#include "Core/Image.hpp"
#include "Core/AssetLoader.hpp"

namespace Panda
{
    class SceneObjectTexture: public BaseSceneObject
    {
        protected:
            std::string m_Name;
            uint32_t m_TexCoordIndex;
            std::shared_ptr<Image> m_pImage;

            std::vector<Matrix4f> m_Transforms;

        public:
            SceneObjectTexture() : BaseSceneObject(SceneObjectType::kSceneObjectTypeTexture), m_TexCoordIndex(0), m_pImage(nullptr) {}
			SceneObjectTexture(const std::string& name) : BaseSceneObject(SceneObjectType::kSceneObjectTypeTexture), m_Name(name), m_TexCoordIndex(0), m_pImage(nullptr) {}
            SceneObjectTexture(uint32_t coord_index, std::shared_ptr<Image>& image) : BaseSceneObject(SceneObjectType::kSceneObjectTypeTexture), m_TexCoordIndex(coord_index), m_pImage(image) {}
            SceneObjectTexture(uint32_t coord_index, std::shared_ptr<Image>&& image) : BaseSceneObject(SceneObjectType::kSceneObjectTypeTexture), m_TexCoordIndex(coord_index), m_pImage(std::move(image)) {}
            SceneObjectTexture(SceneObjectTexture&) = default;
            SceneObjectTexture(SceneObjectTexture&&) = default;
            
            void AddTransform(Matrix4f& matrix) {m_Transforms.push_back(matrix);}
            void SetName(const std::string& name) {m_Name = name;}
            void SetName(std::string&& name) {m_Name = std::move(name);}
            const std::string& GetName() const {return m_Name;}
            void LoadTexture()
            {
                if (!m_pImage)
                {
                    // we should lookup if the texture has been loaded already to prevent
                    // duplicate load. This could be done in Asset Loader Manager.
                    Buffer buf = g_pAssetLoader->SyncOpenAndReadBinary(m_Name.c_str());
                    std::string ext = m_Name.substr(m_Name.find_last_of("."));
                    if (ext == ".jpg" || ext == "jpeg")
                    {
                        JfifParser jfifParser;
                        m_pImage = std::make_shared<Image>(jfifParser.Parse(buf));
                    }
                    else if (ext == ".png")
                    {
                        PngParser pngParser;
                        m_pImage = std::make_shared<Image>(pngParser.Parse(buf));
                    }
                    else if (ext == ".bmp")
                    {
                        BmpParser bmpParser;
                        m_pImage = std::make_shared<Image>(bmpParser.Parse(buf));
                    }
                    else if (ext == ".tga")
                    {
                        TgaParser tgaParser;
                        m_pImage = std::make_shared<Image>(tgaParser.Parse(buf));
                    }
                }
            }

            const Image& GetTextureImage()
            {
                if (!m_pImage)
                {
                    LoadTexture();
                }

                return *m_pImage;
            }

            friend std::ostream& operator<<(std::ostream& out, const SceneObjectTexture& obj);
    };
}