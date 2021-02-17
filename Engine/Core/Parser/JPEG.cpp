#include "JPEG.hpp"
#include "turbojpeg.h"

namespace Panda
{
    size_t JfifParser::ParseScanData(const uint8_t* pScanData, const uint8_t* pDataEnd, Image& img)
    {
        std::vector<uint8_t> scanData;
        size_t scanLength = 0;

        {
            const uint8_t* p = pScanData;

            // scan for scan data buffer size and remove bitstuff
            bool bitstuff = false;
            while (p < pDataEnd && (*p != 0xFF || *(p + 1) == 0x00))
            {
                if (!bitstuff)
                {
                    scanData.push_back(*p);
                }
                else 
                {
                    // ignore it and reset the flag
                    assert(*p == 0x00);
                    bitstuff = false;
                }

                if (*(uint16_t*)p == to_endian_net((uint16_t)0xFF00))
                    bitstuff = true;

                p++;
                scanLength++;
            }

            if (*p == 0xFF && *(p + 1) >= 0xD0 && *(p + 1) <= 0xD7)
            {
                // found restart mark
#if DUMP_DETAILS
                std::cout << "Found RST while scan the ECS." << std::endl;
#endif
            }

#if DUMP_DETAILS
            std::cout << "Size of Scan: " << scanLength << " bytes" << std::endl;
            std::cout << "Size of Scan (after remove bitstuff):" << scanData.size() << " bytes" << std::endl;
            std::cout << "Total MCU count: " << McuCount << std::endl;
#endif
        }

        int16_t previousDC[4]; // 4 is max num of components defined by ITU-T81
        memset(previousDC, 0x00, sizeof(previousDC));

        size_t byteOffset = 0;
        uint8_t bitOffset = 0;

        while(byteOffset < scanData.size() && McuIndex < McuCount)
        {
            #if DUMP_DETAILS
            std::cout << "MCU: " << McuIndex << std::endl;
            #endif
            float block[4][64]; // 4 is max num of components defined by ITU-T81
            memset(&block, 0x00, sizeof(block));

            for (uint8_t i = 0; i < m_ComponentsInFrame; ++i)
            {
                const FRAME_COMPONENT_SPEC_PARAMS& fcsp = m_TableFrameComponentSpec[i];
                #if DUMP_DETAILS
                std::cout << "\tComponent Selector: " << (uint16_t)pScsp[i].ComponentSelector << std::endl;
                std::cout << "\tQuantization Table Destination Selector: " << (uint16_t)fcsp.QuantizationTableDestSelector << std::endl;
                std::cout << "\tDC Entropy Coding table Destination Selector: " << (uint16_t)pScsp[i].DcEntropyCodingTableDestSelector() << std::endl;
                std::cout << "\tAC Entropy Coding table Destination Selector: " << (uint16_t)pScsp[i].AcEntropyCodingTableDestSelector() << std::endl;
                #endif 

                // Decode DC
                uint8_t dcCode = m_TreeHuffman[pScsp[i].DcEntropyCodingTableDestSelector()].DecodeSingleValue(scanData.data(), scanData.size(), &byteOffset, &bitOffset);
                uint8_t dcBitLength = dcCode & 0x0F;
                int16_t dcValue;
                uint32_t tmpValue;

                if (!dcCode)
                {
                    #if DUMP_DETAILS
                    std::cout << "Found EOB when decode DC!" << std::endl;
                    #endif

                    dcValue = 0;
                }
                else
                {
                    if (dcBitLength + bitOffset <= 8)
                    {
                        tmpValue = ((scanData[byteOffset] & ((0x01u << (8 - bitOffset)) - 1)) >> (8 - dcBitLength - bitOffset));
                    }
                    else 
                    {
                        uint8_t bitsInFirstByte = 8 - bitOffset;
                        uint8_t appendFullBytes = (dcBitLength - bitsInFirstByte) / 8;
                        uint8_t bitsInLastByte = dcBitLength - bitsInFirstByte - 8 * appendFullBytes;
                        tmpValue = (scanData[byteOffset] & ((0x01u << (8 - bitOffset)) - 1));
                        for (size_t m = 1; m <= appendFullBytes; ++m)
                        {
                            tmpValue <<= 8;
                            tmpValue += scanData[byteOffset + m];
                        }
                        tmpValue <<= bitsInLastByte;
                        tmpValue += (scanData[byteOffset + appendFullBytes + 1] >> (8 - bitsInLastByte));
                    }

                    // decode dc value;
                    if ((tmpValue >> (dcBitLength - 1)) == 0)
                    {
                        // MSB = 1, turn it to minus value
                        dcValue = -(int16_t)(~tmpValue & ((0x0001u << dcBitLength) - 1));
                    }
                    else 
                    {
                        dcValue = tmpValue;
                    }
                }
                
                // add with previous DC value
                dcValue += previousDC[i];
                // save the value for next DC
                previousDC[i] = dcValue;

                #ifdef DUMP_DETAILS
                printf("DC Code: %x\n", dcCode);
                printf("DC Bit Length:%d\n", dcBitLength);
                printf("DC Value: %d\n", dcValue);
                #endif

                block[i][0] = dcValue;

                // forward pointers to end of DC
                bitOffset += dcBitLength;
                while(bitOffset >= 8)
                {
                    bitOffset -= 8;
                    byteOffset++;
                }

                // Decode AC
                int32_t acIndex = 1;
                while(byteOffset < scanData.size() && acIndex < 64)
                {
                    uint8_t acCode = m_TreeHuffman[2 + pScsp[i].AcEntropyCodingTableDestSelector()].DecodeSingleValue(scanData.data(), scanData.size(), &byteOffset, &bitOffset);

                    if (!acCode)
                    {
                        #if DUMP_DETAILS
                        std::cout << "Found EOB when decode AC!" << std::endl;
                        #endif
                        break;
                    }
                    else if (acCode == 0xF0)
                    {
                        #if DUMP_DETAILS
                        std::cout << "Found ZRL when decode AC!" << std::endl;
                        #endif
                        acIndex += 16;
                        continue;
                    }

                    uint8_t acZeroLength = acCode >> 4;
                    acIndex += acZeroLength;
                    uint8_t acBitLength = acCode & 0x0F;
                    int16_t acValue;

                    if (acBitLength + bitOffset <= 8)
                    {
                        tmpValue = ((scanData[byteOffset] & ((0x01u << (8 - bitOffset)) - 1)) >> (8 - acBitLength - bitOffset));
                    }
                    else 
                    {
                        uint8_t bitsInFirstByte = 8 - bitOffset;
                        uint8_t appendFullBytes = (acBitLength - bitsInFirstByte) / 8;
                        uint8_t bitsInLastByte = acBitLength - bitsInFirstByte - 8 * appendFullBytes;
                        tmpValue = (scanData[byteOffset] & ((0x01u << (8 - bitOffset)) - 1));
                        for (size_t m = 1; m <= appendFullBytes; ++m)
                        {
                            tmpValue <<= 8;
                            tmpValue += scanData[byteOffset + m];
                        }
                        tmpValue <<= bitsInLastByte;
                        tmpValue += (scanData[byteOffset + appendFullBytes + 1] >> (8 - bitsInLastByte));
                    }

                    // decode ac value
                    if ((tmpValue >> (acBitLength - 1)) == 0)
                    {
                        // MSB = 1, turn it to minus value
                        acValue = -(int16_t)(~tmpValue & ((0x0001u << acBitLength) - 1));
                    }
                    else
                    {
                        acValue = tmpValue;
                    }

                    #ifdef DUMP_DETAILS
                    printf("AC Code: %x\n", acCode);
                    printf("AC Bit Length: %d\n", acBitLength);
                    printf("AC Value: %d\n", acValue);
                    #endif

                    int32_t index = m_ZigzagIndex[acIndex];
                    block[i][(index >> 3) * 8 + (index & 0x07)] = acValue;

                    // forward pointers to end of AC
                    bitOffset += acBitLength;
                    while(bitOffset >= 8)
                    {
                        bitOffset -= 8;
                        byteOffset++;
                    }

                    acIndex++;
                }

                #ifdef DUMP_DETAILS
                printf("Extracted Component[%d] 8x8 block: \n", i);
                for (size_t _i = 0; _i < 64; ++_i)
                {

                    if (_i != 0 && (_i % 8 == 0))
                        std::cout << std::endl;
                    std::cout << block[i][_i] << ",";
                }
                std::cout << std::endl << std::endl;
                #endif

                for (size_t _i = 0; _i < 64; ++_i)
                {
                    block[i][_i] *= m_TableQuantization[fcsp.QuantizationTableDestSelector][_i];
                }
                
#ifdef DUMP_DETAILS
                std::cout << "After Quantization: " << std::endl;
                for (size_t _i = 0; _i < 64; ++_i)
                {

                    if (_i != 0 && (_i % 8 == 0))
                        std::cout << std::endl;
                    std::cout << block[i][_i] << ",";
                }
                std::cout << std::endl << std::endl;
#endif

                block[i][0] += 1024.0f; // level shift. same as +128 to each element after IDCT
                IDCT8x8(block[i], block[i]);
                #ifdef DUMP_DETAILS
                std::cout << "After IDCT: " << std::endl;
                for (size_t _i = 0; _i < 64; ++_i)
                {

                    if (_i != 0 && (_i % 8 == 0))
                        std::cout << std::endl;
                    std::cout << block[i][_i] << ",";
                }
                std::cout << std::endl << std::endl;
                #endif
            }

            assert(m_ComponentsInFrame <= 4);

            YCbCrf ycbcr;
            RGBf rgb;
            int mcuIndexX = McuIndex % McuCountX;
            int mcuIndexY = McuIndex / McuCountX;
            uint8_t* pBuf;

            for (size_t i = 0; i < 8; ++i)
                for (size_t j = 0; j < 8; ++j)
                {
                    for (size_t k = 0; k < m_ComponentsInFrame; ++k)
                    {
                        ycbcr.data[k] = block[k][i * 8 + j];
                    }

                    pBuf = reinterpret_cast<uint8_t*>(img.Data) + (img.Pitch * (mcuIndexY * 8 + i) + (mcuIndexX * 8 + j) * (img.BitCount >> 3));
                    rgb = ConvertYCbCr2RGB(ycbcr);
                    reinterpret_cast<R8G8B8A8Unorm*>(pBuf)->data[0] = (uint8_t)rgb[0];
                    reinterpret_cast<R8G8B8A8Unorm*>(pBuf)->data[1] = (uint8_t)rgb[1];
                    reinterpret_cast<R8G8B8A8Unorm*>(pBuf)->data[2] = (uint8_t)rgb[2];
                    reinterpret_cast<R8G8B8A8Unorm*>(pBuf)->data[3] = 255;
                }

            McuIndex++;

            if(m_RestartInterval != 0 && (McuIndex % m_RestartInterval == 0))
            {
                if (bitOffset)
                {
                    // finish current byte
                    bitOffset = 0;
                    byteOffset++;
                }
                assert(byteOffset == scanData.size());
                break;
            }
        }
        
        return scanLength;
    }

	Image JfifParser::Parse(Buffer& buf)
	{
		Image img;
		const char *subsampName[TJ_NUMSAMP] = {
			"4:4:4", "4:2:2", "4:2:0", "Grayscale", "4:4:0", "4:1:1"
		};

		const char *colorspaceName[TJ_NUMCS] = {
			"RGB", "YCbCr", "GRAY", "CMYK", "YCCK"
		};

		// Init decompress instance
		tjhandle tjInstance = tjInitDecompress();
		assert(tjInstance != NULL);

		int32_t width, height;
		int32_t inSubsamp, inColorspace;
		if (tjDecompressHeader3(tjInstance, buf.GetData(), buf.GetDataSize(), &width, &height,
			&inSubsamp, &inColorspace) < 0)
		{
			assert(false);
			printf("Read jpeg file header failed. Please check the file!");
			return img;
		}

		printf("Image:  %d x %d pixels, %s subsampling, %s colorspace\n",
			width, height, subsampName[inSubsamp], colorspaceName[inColorspace]);

		int32_t pixelFormat = TJPF_RGBA;
		img.Width = width;
		img.Height = height;
		img.BitCount = 32;
		img.Pitch = img.Width * (img.BitCount >> 3);
		img.DataSize = width * height * tjPixelSize[pixelFormat];
		img.Data = g_pMemoryManager->Allocate(img.DataSize);
		if (tjDecompress2(tjInstance, buf.GetData(), buf.GetDataSize(), (unsigned char*)img.Data, img.Width,
			0, img.Height, pixelFormat, 0) < 0)
		{
			assert(false);
			printf("Decompress jpeg file failed. Please check the file!");
			return img;
		}

		tjDestroy(tjInstance);

		return img;
	}

    Image JfifParser::Parser1(Buffer& buf)
    {
        Image img;

        const uint8_t* pData = buf.GetData();
        const uint8_t* pDataEnd = buf.GetData() + buf.GetDataSize();

        const JFIF_FILEHEADER* pFileHeader = reinterpret_cast<const JFIF_FILEHEADER*>(pData);
        pData += sizeof(JFIF_FILEHEADER);
        if (pFileHeader->SOI == to_endian_net((uint16_t)0xFFD8))
        {
            std::cout << "Asset if JPEG file" << std::endl;

            while(pData < pDataEnd)
            {
                size_t scanLength = 0;

                const JPEG_SEGMENT_HEADER* pSegmentHeader = reinterpret_cast<const JPEG_SEGMENT_HEADER*>(pData);

                switch(to_endian_native(pSegmentHeader->Marker))
                {
                    case 0xFFC0:
                    case 0xFFC2:
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        if (to_endian_native(pSegmentHeader->Marker) == 0xFFC0)
                            std::cout << "Start of Frame0 (baseline DCT)" << std::endl;
                        else 
                            std::cout << "Start of Frame2 (progressive DCT)" << std::endl;

                        std::cout << "-----------------------" << std::endl;

                        const FRAME_HEADER* pFrameHeader = reinterpret_cast<const FRAME_HEADER*>(pData);
                        m_SamplePrecision = pFrameHeader->SamplePrecision;
                        m_Lines = to_endian_native(pFrameHeader->NumOfLines);
                        m_SamplesPerLine = to_endian_native(pFrameHeader->NumOfSamplesPerLine);
                        m_ComponentsInFrame = pFrameHeader->NumOfComponentsInFrame;
                        McuIndex = 0;
                        McuCountX = ((m_SamplesPerLine + 7) >> 3);
                        McuCountY = ((m_Lines + 7) >> 3);
                        McuCount = McuCountX * McuCountY;

                        std::cout << "Sample Precision: " << m_SamplePrecision << std::endl;
                        std::cout << "Num of Lines: " << m_Lines << std::endl;
                        std::cout << "Num of Samples per line: " << m_SamplesPerLine << std::endl;
                        std::cout << "Num of Components In Frame: " << m_ComponentsInFrame << std::endl;
                        std::cout << "Totale MCU count: " << McuCount << std::endl;

                        const uint8_t* pTmp = pData + sizeof(FRAME_HEADER);
                        const FRAME_COMPONENT_SPEC_PARAMS* pFcsp = reinterpret_cast<const FRAME_COMPONENT_SPEC_PARAMS*>(pTmp);
                        for (uint8_t i = 0; i < pFrameHeader->NumOfComponentsInFrame; ++i)
                        {
                            std::cout << "\tComponent Identifier: " << (uint16_t)pFcsp->ComponentIdentifier << std::endl;
                            std::cout << "\tHorizontal Sampling Factor: " << (uint16_t)pFcsp->HorizontalSamplingFactor() << std::endl;
                            std::cout << "\tVertical Sampling Factor:" << (uint16_t)pFcsp->VerticalSamplingFactor() << std::endl;
                            std::cout << "\tQuantization Table Destination Selector: " << (uint16_t)pFcsp->QuantizationTableDestSelector << std::endl;
                            std::cout << std::endl;
                            m_TableFrameComponentSpec.push_back(*pFcsp);
                            pFcsp++;
                        }

                        img.Width = m_SamplesPerLine;
                        img.Height = m_Lines;
                        img.BitCount = 32;
                        img.Pitch = McuCountX * 8 * (img.BitCount >> 3);
                        img.DataSize = img.Pitch * McuCountY * 8; //* (img.BitCount >> 3);
                        img.Data = g_pMemoryManager->Allocate(img.DataSize);

                        std::cout << std::endl;
                        pData += to_endian_native(pSegmentHeader->Length) + 2; // length of marker
                    }
                    break;
                    case 0xFFC4: // Define Huffman Table(s)
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        std::cout << "Define Huffman Table(s)" << std::endl;
                        std::cout << "-----------------------" << std::endl;

                        size_t segmentLength = to_endian_native(pSegmentHeader->Length) - 2;

                        const uint8_t* pTmp = pData + sizeof(JPEG_SEGMENT_HEADER);

                        while (segmentLength > 0)
                        {
                            const HUFFMAN_TABLE_SPEC* pHtable = reinterpret_cast<const HUFFMAN_TABLE_SPEC*>(pTmp);
                            std::cout << "Table Class: " << pHtable->TableClass() << std::endl;
                            std::cout << "Destination Identifier: " << pHtable->DestinationIdentifier() << std::endl;

                            const uint8_t* pCodeValueStart = reinterpret_cast<const uint8_t*>(pHtable) + sizeof(HUFFMAN_TABLE_SPEC);
                            auto numSymbo = m_TreeHuffman[(pHtable->TableClass() << 1) | pHtable->DestinationIdentifier()].PopulateWithHuffmanTable(pHtable->NumOfHuffmanCodes, pCodeValueStart);

                            #ifdef DUMP_DETAILS
                            m_TreeHuffman[(pHtable->TableClass() << 1) | pHtable->DestinationIdentifier()].Dump();
                            #endif

                            size_t processedLength = sizeof(HUFFMAN_TABLE_SPEC) + numSymbo;
                            pTmp += processedLength;
                            segmentLength -= processedLength;
                        }
                        std::cout << std::endl;
                        pData += to_endian_native(pSegmentHeader->Length) + 2; // length of the marker
                    }
                    break;
                    case 0xFFDB: // Define Quantization Table(s)
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        std::cout << "Define Quantization Table(s) " << std::endl;
                        std::cout << "-----------------------------" << std::endl;

                        size_t segmentLength = to_endian_native(pSegmentHeader->Length) - 2;
                        const uint8_t* pTmp = pData + sizeof(JPEG_SEGMENT_HEADER);

                        while (segmentLength > 0)
                        {
                            const QUANTIZATION_TABLE_SPEC* pQtable = reinterpret_cast<const QUANTIZATION_TABLE_SPEC*>(pTmp);
                            std::cout << "Element Precision: " << pQtable->ElementPrecision() << std::endl;
                            std::cout << "Destination Identifier: " << pQtable->DestinationIdentifier() << std::endl;

                            const uint8_t* pElementDataStart = reinterpret_cast<const uint8_t*>(pQtable) + sizeof(QUANTIZATION_TABLE_SPEC);

                            for (size_t i = 0; i < 64; ++i)
                            {
                                int32_t index = m_ZigzagIndex[i];
                                if (pQtable->ElementPrecision() == 0)
                                    m_TableQuantization[pQtable->DestinationIdentifier()][(index >> 3) * 8 + (index & 0x7)] = pElementDataStart[i];
                                else 
                                    m_TableQuantization[pQtable->DestinationIdentifier()][(index >> 3) * 8 + (index & 0x7)] = to_endian_native(*((uint16_t*)pElementDataStart + i));
                            }
                            #ifdef DUMP_DETAILS
                            std::cout << std::endl;
                            for (size_t i = 0; i < 64; ++i)
                            {
                                if (i != 0 && (i % 8 == 0))
                                    std::cout << std::endl;
                                std::cout << m_TableQuantization[pQtable->DestinationIdentifier()][i] << ",";
                            }
                            std::cout << std::endl << std::endl;
                            #endif
                            size_t processedLength = sizeof(QUANTIZATION_TABLE_SPEC) + 64 * (pQtable->ElementPrecision() + 1);
                            pTmp += processedLength;
                            segmentLength -= processedLength;
                        }
                        std::cout << std::endl;
                        pData += to_endian_native(pSegmentHeader->Length) + 2; // lengh of the marker
                    }
                    break;
                    case 0xFFDD:
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        std::cout << "Define Restart Interval" << std::endl;
                        std::cout << "-----------------------" << std::endl;

                        RESTART_INTERVAL_DEF* pRestartHeader = (RESTART_INTERVAL_DEF*)pData;
                        m_RestartInterval = to_endian_native((uint16_t)pRestartHeader->RestertInterval);
                        std::cout << "Restart interval: " << m_RestartInterval << std::endl;
                        std::cout << std::endl;
                        pData += to_endian_native(pSegmentHeader->Length) + 2; // length of the marker
                    }
                    break;
                    case 0xFFDA:
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        std::cout << "Start Of Scan" << std::endl;
                        std::cout << "----------------------" << std::endl;

                        SCAN_HEADER* pScanHeader = (SCAN_HEADER*)pData;
                        std::cout << "Image Components in Scan: " << (uint16_t)pScanHeader->NumOfComponents << std::endl;
                        assert(pScanHeader->NumOfComponents == m_ComponentsInFrame);

                        const uint8_t* pTmp = pData + sizeof(SCAN_HEADER);
                        pScsp = reinterpret_cast<const SCAN_COMPONENT_SPEC_PARAMS*>(pTmp);

                        const uint8_t* pScanData = pData + to_endian_native((uint16_t)pScanHeader->Length) + 2;

                        scanLength = ParseScanData(pScanData, pDataEnd, img);
                        std::cout << std::endl;
                        pData += to_endian_native(pSegmentHeader->Length) + 2 + scanLength; // length of the marker
                    }
                    break;
                    case 0xFFD0:
                    case 0xFFD1:
                    case 0xFFD2:
                    case 0xFFD3:
                    case 0xFFD4:
                    case 0xFFD5:
                    case 0xFFD6:
                    case 0xFFD7:
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        #if DUMP_DETAILS
                        std::cout << "Restart of Scan" << std::endl;
                        std::cout << "----------------------" << std::endl;
                        #endif

                        const uint8_t* pScanData = pData + 2;
                        scanLength = ParseScanData(pScanData, pDataEnd, img);
                        std::cout << std::endl;
                        pData += 2 + scanLength; // length of the marker
                    }
                    break;
                    case 0xFFD9:
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "End of Scan" << std::endl;
                        std::cout << "------------------------" << std::endl;
                        std::cout << std::endl;
                        pData += 2; // length of marker
                    }
                    break;
                    case 0xFFE0: // Application specific
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        const APP0* pApp0 = reinterpret_cast<const APP0*>(pData);
                        switch(to_endian_native(*(uint32_t*)pApp0->Identifier))
                        {
                            case "JFIF\0"_u32:
                            {
                                const JFIF_APPO* pJfifApp0 = reinterpret_cast<const JFIF_APPO*>(pApp0);
                                std::cout << "JFIF-APP0" << std::endl;
                                std::cout << "--------------------------" << std::endl;
                                std::cout << "JFIF Version: " << (uint16_t)pJfifApp0->MajorVersion << "." 
                                    << (uint16_t)pJfifApp0->MinorVersion << std::endl;
                                std::cout << "Density Units: " <<
                                    ((pJfifApp0->DensityUnits == 0)? "No units" : 
                                        ((pJfifApp0->DensityUnits == 1)? "Pixel per inch" : "Pixel per centimeter"))
                                    << std::endl;
                                std::cout << "Density: " << to_endian_native(pJfifApp0->Xdensity) << "*"
                                    << to_endian_native(pJfifApp0->Ydensity) << std::endl;
                                if (pJfifApp0->Xthumbnail && pJfifApp0->Ythumbnail)
                                {
                                    std::cout << "Thumbnail Dimesions[w * h]: " << (uint16_t)pJfifApp0->Xthumbnail << "*"
                                        << (uint16_t)pJfifApp0->Ythumbnail << std::endl;
                                }
                                else 
                                {
                                    std::cout << "No thumbnail defined in JFIF-APP0 segment!" << std::endl;
                                }
                                
                            }
                            break;
                            case "JFXX\0"_u32:
                            {
                                const JFXX_APP0* pJfxxApp0 = reinterpret_cast<const JFXX_APP0*>(pApp0);
                                std::cout << "Thumbnail Format: ";
                                switch(pJfxxApp0->ThumbnailFormat)
                                {
                                    case 0x10:
                                        std::cout << "JPEG format";
                                        break;
                                    case 0x11:
                                        std::cout << "1 byte per pixel palettized format";
                                        break;
                                    case 0x13:
                                        std::cout << "3 bytes per pixel RGB format";
                                        break;
                                    default:
                                        std::printf("Unrecognized Thumbnail Format: %x\n", pJfxxApp0->ThumbnailFormat);
                                        break;
                                }
                            }
                            break;
                            default:
                                std::cout << "Ignore Unrecognized APP0 segment" << std::endl;
                        }
                        std::cout << std::endl;
                        pData += to_endian_native(pSegmentHeader->Length) + 2;
                    }
                    break;
                    case 0xFFFE:
                    {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                        std::cout << "Text Comment" << std::endl;
                        std::cout << "---------------------" << std::endl;
                        std::cout << std::endl;
                        pData += to_endian_native(pSegmentHeader->Length) + 2;
                    }
                    break;
                    default:
                        {
#if DUMP_DETAILS
                        std::cout << "======================" << std::endl;
#endif
                        std::cout << "Segment Length: " << to_endian_native(pSegmentHeader->Length) << " bytes" << std::endl;
                            std::printf("Ignore unrecognized Segment. Marker = %0x\n\n", to_endian_native(pSegmentHeader->Marker));
                            std::cout << std::endl;
                            pData += to_endian_native(pSegmentHeader->Length) + 2; // length of the marker
                        }
                        break;
                }
            }
        }
        else 
            std::cout << "File is not a JPEG file!" << std::endl;

        return img;
    }   
}