#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>
#include <string>
#include "BaseSceneNode.hpp"

namespace Panda
{    
    std::ostream& operator<<(std::ostream& out, const BaseSceneNode& node)
    {
		static thread_local int32_t indent = 0;
		indent++;

		out << std::string(indent, ' ') << "Scene Node" << std::endl;
		out << std::string(indent, ' ') << "----------" << std::endl;
		out << std::string(indent, ' ') << "Name: " << node.m_Name << std::endl;
		node.Dump(out);
		out << std::endl;

		for (auto subNode : node.m_Children)
		{
			out << *subNode << std::endl;
		}

		for (auto subTransform : node.m_Transforms)
		{
			out << *subTransform << std::endl;
		}

		indent--;

        return out;
    }
}
