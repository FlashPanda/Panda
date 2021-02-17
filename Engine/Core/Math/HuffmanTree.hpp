#pragma once

#include <memory>
#include <vector>
#include "Core/portable.hpp"
#include "MathUtility.hpp"

namespace Panda
{
    template<typename T>
    class HuffmanNode 
    {
        protected:
            T m_Value;
            std::shared_ptr<HuffmanNode<T>> m_pLeft;
            std::shared_ptr<HuffmanNode<T>> m_pRight;
            bool m_IsLeaf = false;

        public:
            HuffmanNode() = default;
            HuffmanNode(T value) : m_Value(value), m_IsLeaf(true) {}
            ~HuffmanNode() = default;

            HuffmanNode(HuffmanNode&) = default;
            HuffmanNode(HuffmanNode&&) = default;
            HuffmanNode& operator=(HuffmanNode&) = default;
            HuffmanNode& operator=(HuffmanNode&&) = default;

            bool IsLeaf() const {return m_IsLeaf;}
            void SetLeft(std::shared_ptr<HuffmanNode> pNode) {m_pLeft = pNode;}
            void SetRight(std::shared_ptr<HuffmanNode> pNode) {m_pRight = pNode;}
            const std::shared_ptr<HuffmanNode<T>> GetLeft() const {return m_pLeft;}
            const std::shared_ptr<HuffmanNode<T>> GetRight() const {return m_pRight;}
            void SetValue(T value) {m_Value = value; m_IsLeaf = true;}
            T GetValue() const {return m_Value;}
    };

    template<typename T>
    class HuffmanTree
    {
        protected:
            std::shared_ptr<HuffmanNode<T>> m_pRoot;

        private:
            void RecursiveDump(const std::shared_ptr<HuffmanNode<T>>& pNode, std::string bitStream)
            {
                if (pNode)
                {
                    if (pNode->IsLeaf())
                    {
                        printf("%20s | %x\n", bitStream.c_str(), pNode->GetValue());
                    }
                    else
                    {
                        RecursiveDump(pNode->GetLeft(), bitStream + "0");
                        RecursiveDump(pNode->GetRight(), bitStream + "1");
                    }
                    
                }
            }

        public:
            HuffmanTree() { m_pRoot = std::make_shared<HuffmanNode<uint8_t>>();}

            size_t PopulateWithHuffmanTable(const uint8_t numOfCodes[16], const uint8_t* codeValues)
            {
                int32_t numSymbo = 0;
                for (size_t i = 0; i < 16; ++i)
                {
                    numSymbo += numOfCodes[i];
                }
                const uint8_t* pCodeValueEnd = codeValues + numSymbo - 1;

                std::queue<std::shared_ptr<HuffmanNode<uint8_t>>> nodeQueue;
                bool foundBottom = false;

                for (int32_t i = 15; i >= 0; --i)
                {
                    int32_t l = numOfCodes[i];
                    if (!foundBottom)
                    {
                        if (l == 0)
                        {
                            // simply move to upper layer
                            continue;
                        }
                        else
                        {
                            foundBottom = true;
                        }
                    }

                    auto childrenCount = nodeQueue.size();

                    if (l)
                    {
                        // create leaf node for level i
                        pCodeValueEnd = pCodeValueEnd - l;
                        const uint8_t* pCodeValue = pCodeValueEnd + 1;

                        for(int j = 0; j < l; ++j)
                        {
                            auto pNode = std::make_shared<HuffmanNode<uint8_t>>(*pCodeValue++);
                            nodeQueue.push(pNode);
                        }
                    }

                    // create non-leaf node and append children
                    while(childrenCount > 0)
                    {
                        auto pNode = std::make_shared<HuffmanNode<uint8_t>>();
                        auto pLeftNode = nodeQueue.front();
                        nodeQueue.pop();
                        pNode->SetLeft(pLeftNode);
                        childrenCount--;

                        if (childrenCount > 0)
                        {
                            auto pRightNode = nodeQueue.front();
                            nodeQueue.pop();
                            pNode->SetRight(pRightNode);
                            childrenCount--;
                        }

                        nodeQueue.push(pNode);
                    }
                }

                // now append to the root node
                assert(nodeQueue.size() <= 2 && nodeQueue.size() > 0);
                auto pLeftNode = nodeQueue.front();
                nodeQueue.pop();
                m_pRoot->SetLeft(pLeftNode);
                if (!nodeQueue.empty())
                {
                    auto pRightNode = nodeQueue.front();
                    nodeQueue.pop();
                    m_pRoot->SetRight(pRightNode);
                }

                return numSymbo;
            }

            T DecodeSingleValue(const uint8_t* encodedStream, const size_t encodeStreamLength, 
                size_t* byteOffset, uint8_t* bitOffset)
            {
                T res = 0;
                std::shared_ptr<HuffmanNode<T>> pNode = m_pRoot;
                for (size_t i = *byteOffset; i < encodeStreamLength; ++i)
                {
                    uint8_t data = encodedStream[i];
                    for (size_t j = *bitOffset; j < 8; ++j)
                    {
                        uint8_t bit = (data & (0x1 << (7 - j))) >> (7 - j);
                        if (bit == 0)
                            pNode = pNode->GetLeft();
                        else
                            pNode = pNode->GetRight();

                        assert(pNode);

                        if (pNode->IsLeaf())
                        {
                            // move pointers to next bit
                            if (j == 7)
                            {
                                *bitOffset = 0;
                                *byteOffset = i + 1;
                            }
                            else
                            {
                                *bitOffset = j + 1;
                                *byteOffset = i;
                            }
                            res = pNode->GetValue();
                            return res;
                        }
                    }
                    *bitOffset = 0;
                }

                // decode failed
                *byteOffset = -1;
                *bitOffset = -1;

                return res;
            }

            std::vector<T> Decode(const uint8_t* encodedStream, const size_t encodedStreamLength)
            {
                std::vector<T> res;
                std::shared_ptr<HuffmanNode<T>> pNode = m_pRoot;

                for (size_t i = 0; i < encodedStreamLength; ++i)
                {
                    uint8_t data = encodedStream[i];
                    for (int j = 0; j < 8; ++j)
                    {
                        uint8_t bit = (data & (0x1 << (7 - j))) >> (7 - j);
                        if (bit == 0)
                            pNode = pNode->GetLeft();
                        else
                            pNode = pNode->GetRight();

                        assert(pNode);

                        if (pNode->IsLeaf())
                        {
                            res.push_back(pNode->GetValue());
                            pNode = m_pRoot;
                        }
                    }
                }

                return res;
            }

            void Dump()
            {
                std::string bitStream;
				RecursiveDump(m_pRoot, bitStream);
                std::cout << std::endl;
            }
    };
}