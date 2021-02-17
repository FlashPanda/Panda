#pragma once

/**
 * wrap the sync read file functions for asyn load files.
 */ 

#include <cstdio>
#include <string>
#include <utility>
#include <vector>
#include "Interface/IRuntimeModule.hpp"
#include "Buffer.hpp"

namespace Panda
{
    class AssetLoader : public IRuntimeModule
    {
        public:
            virtual ~AssetLoader() {}

            virtual int Initialize();
            virtual void Finalize();

            virtual void Tick();

            typedef void* AssetFilePtr;

            enum AssetOpenMode
            {
                PANDA_OPEN_TEXT     = 0,    // open in text mode
                PANDA_OPEN_BINARY   = 1    // open in binary mode
            };

            enum AssetSeekBase
            {
                PANDA_SEEK_SET = 0, // seek set
                PANDA_SEEK_CUR = 1, // seek cur
                PANDA_SEEK_END = 2  // seek end
            };

            bool AddSearchPath(const char* path);
            bool RemoveSearchPath(const char* path);

            virtual bool FileExists(const char* filePath);

            virtual AssetFilePtr OpenFile(const char* name, AssetOpenMode mode);

            virtual Buffer SyncOpenAndReadText(const char* filePath);

			virtual Buffer SyncOpenAndReadBinary(const char* filePath);

            virtual size_t SyncRead(const AssetFilePtr& fp, Buffer& buf);

            virtual void CloseFile(AssetFilePtr& fp);

            virtual size_t GetSize(const AssetFilePtr& fp);

            virtual int32_t Seek(AssetFilePtr fp, long offset, AssetSeekBase where);

            inline std::string SyncOpenAndReadFileToString (const char* fileName)
            {
                std::string result;
                Buffer buffer = SyncOpenAndReadText(fileName);
                if (buffer.GetDataSize())
                {
                    char* content = reinterpret_cast<char*>(buffer.GetData());
                    if (content)
                    {
                        result = std::string(std::move(content)); // move content and clear content values, since C++11
                    }
                }

                return result;
            }

        private:
            std::vector<std::string> m_SearchPath;
    };

    extern AssetLoader* g_pAssetLoader;
}
