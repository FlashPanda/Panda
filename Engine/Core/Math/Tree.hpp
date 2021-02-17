#pragma once
#include <iostream>
#include <vector>
#include "MathUtility.hpp"

namespace Panda
{
    class TreeNode
    {
        protected:
            TreeNode* m_Parent;
            std::vector<std::shared_ptr<TreeNode>> m_Children;

        protected:
            virtual void Dump(std::ostream& out) const {}

        public:
            virtual ~TreeNode() {}

			virtual TreeNode* GetParent()
			{
				return m_Parent;
			}
            virtual void AppendChild(std::shared_ptr<TreeNode>&& sub_node)
            {
                sub_node->m_Parent = this;
                m_Children.push_back(std::move(sub_node));
            }

            friend std::ostream& operator<<(std::ostream& out, const TreeNode& node)
            {
                static thread_local int32_t indent = 0;
                indent ++;

                out << std::string(indent, ' ') << "Tree Node" << std::endl;
                out << std::string(indent, ' ') << "--------------" << std::endl;
                node.Dump(out);
                out << std::endl;

                for (const std::shared_ptr<TreeNode>& sub_node : node.m_Children)
                {
                    out << *sub_node << std::endl;
                }

                indent--;

                return out;
            }
    };

    template <typename T>
    class BinaryTreeNode
    {
    public:
        BinaryTreeNode() :m_Parent(nullptr), m_Left(nullptr), m_Right(nullptr), m_Key(T(0)) {}
		BinaryTreeNode(T key) : m_Parent(nullptr), m_Left(nullptr), m_Right(nullptr), m_Key(key) {}
        BinaryTreeNode(BinaryTreeNode* p, BinaryTreeNode* l, BinaryTreeNode* r, T key)
            : m_Parent(p), m_Left(l), m_Right(r), m_Key(key) {}
        virtual ~BinaryTreeNode() {}

    public:
        BinaryTreeNode* Parent()
        {
            return m_Parent;
        }
        void Parent(BinaryTreeNode* p)
        {
            m_Parent = p;
        }

        BinaryTreeNode* Left()
        {
            return m_Left;
        }
        void Left(BinaryTreeNode* l)
        {
            m_Left = l;
        }

        BinaryTreeNode* Right()
        {
            return m_Right;
        }
        void Right(BinaryTreeNode* r)
        {
            m_Right = r;
        }

        T Key()
        {
            return m_Key;
        }
        void Key(T key)
        {
            m_Key = key;
        }

    protected:
        virtual void Dump(std::ostream& out) const {}

    protected:
        BinaryTreeNode* m_Parent;
        BinaryTreeNode* m_Left;
        BinaryTreeNode* m_Right;
        T m_Key;
    };

    template <typename T>
    class BinarySearchTree
    {
    public:
        BinarySearchTree() : m_Root(nullptr) {}
        BinarySearchTree(BinaryTreeNode<T>* root) : m_Root(root) {}
        virtual ~BinarySearchTree() {}

    public:
        bool IsEmpty() {return m_Root == nullptr?;}
		void Free()
		{
			Free(m_Root);
		}
		BinaryTreeNode<T>* Root()
		{
			return m_Root;
		}
        void Insert(BinaryTreeNode<T>* v)
        {
            // find parent
            BinaryTreeNode<T>* parent = nullptr;
            BinaryTreeNode<T>* p = m_Root;
            while (p != nullptr)
            {
                parent = p;
                if (v->Key() <= p->Key())
                    p = p->Left();
                else 
                    p = p->Right();
            }

            // insert it
            v->Parent(parent);
            if(parent == nullptr) // The tree is empty
                m_Root = v;
            else if (v->Key() <= parent->Key())
                parent->Left(v);
            else 
                parent->Right(v);
        }

        BinaryTreeNode<T>* Search(BinaryTreeNode<T>* p, T key)
        {
            if (p == nullptr || key == p->Key())
                return p;
            else if (key <= p->Key())
                return Search(p->Left(), key);
            else
                return Search(p->Right(), key);
        }

        BinaryTreeNode<T>* Minimum(BinaryTreeNode<T>* p)
        {
            BinaryTreeNode<T>* t = p;
            while (t->Left() != nullptr)
                t = t->Left();
            return t;
        }

        BinaryTreeNode<T>* Maximum(BinaryTreeNode<T>* p)
        {
            BinaryTreeNode<T>* t = p;
            while (t->Right() != nullptr)
                t = t->Right();
            return t;
        }

        BinaryTreeNode<T>* Successor(BinaryTreeNode<T>* p)
        {
            // If the node has a right child, the successor should be 
            // the minimum of the right sub-tree
            if (p->Right() != nullptr)
                return Minimum(p->Right());

            // If the node doesn't have right child, the successor should be
            // the first node treated its father as left child, something like that
            // Or, we will get a nullptr.
            BinaryTreeNode<T>* t1 = p;
            BinaryTreeNode<T>* t2 = p->Parent();
            while (t2 != nullptr && t2->Right() == t1)
            {
                t1 = t2;
                t2 = t2->Parent();
            }

            return t2;
        }

        BinaryTreeNode<T>* Predecessor(BinaryTreeNode<T>* p)
        {
            // If the node has a left child, the successor should be 
            // the maximum of the left sub-tree
            if (p->Left() != nullptr)
                return Maximum(p->Left());

            // If the node doesn't have a left child, the predecessor should be
            // the first node treated by its father as right child, something like that
            // Or, we will get a nullptr.
            BinaryTreeNode<T>* t1 = p;
            BinaryTreeNode<T>* t2 = p->Parent();
            while (t2 != nullptr && t2->Left() == t1)
            {
                t1 = t2;
                t2 = t2->Parent();
            }

            return t2;
        }

        void Transplant(BinaryTreeNode<T>* src, BinaryTreeNode<T>* dst)
        {
            if (src->Parent() == nullptr)
                m_Root = dst;
            else if (src == src->Parent()->Left())
                src->Parent()->Left(dst);
            else 
                src->Parent()->Right(dst);
            if (dst != nullptr)
                dst->Parent(src->Parent());
        }

        void Delete(BinaryTreeNode<T>* src)
        {
            if (src->Left() == nullptr)
                Transplant(src, src->Right()); // It doesn't matter whether right child is or is not nullptr.
            else if (src->Right() == nullptr)
                Transplant(src, src->Left());
            else
            {
                // We need to find a successor to replace _src_
                BinaryTreeNode<T>* s = Minimum(src->Right());
                if (s->Parent() != src)
                {
                    // Use the right child to replace _s_
                    // Cause _s_ might contains a right child.
                    Transplant(s, s->Right());

                    // fill the right child of _s_ with the right child of _src_
                    s->Right(src->Right());
                    s->Right()->Parent(s);
                }
                Transplant(src, s);
                s->Left(src->Left());
                s->Left()->Parent(s);
            }
        }

		void DumpToInorderList(std::vector<T>& list)
		{
			InorderWalk(list, m_Root);
		}

        void DumpWithInorderWalk(std::ostream& out)
        {
            out << "【Inorder Walk】 The tree is: \n";
            InorderWalk(out, m_Root);
			out << std::endl << std::endl;
        }

		void DumpWithPreorderWalk(std::ostream& out)
		{
			out << "【Preorder Walk】The tree is : \n";
			PreorderWalk(out, m_Root);
			out << std::endl << std::endl;
		}

		void DumpStructureToFile(std::fstream& out)
		{
			int32_t h = Height();
			assert(h != 0);
			if (h == 0)
				return;

			std::vector<std::vector<int32_t>> levels;
			std::vector<std::vector<int32_t>> indices;
			std::vector<BinaryTreeNode<int32_t>*> currentLevel;
			currentLevel.push_back(m_Root);
			for (int32_t i = 1; i <= h; ++i)
			{
				// The indices of current level
				std::vector<int32_t> ci;
				std::vector<int32_t> li;
				if (indices.size() > 0)
					li = indices[indices.size() - 1];
				if (li.size() == 0) // We are now dealing with root level.
					ci.push_back(std::pow(2, h - 1) - 1);
				else
				{
					int32_t sbc = std::pow(2, h - i);// side bar count of last level 
					for (int32_t j = 0; j < li.size(); ++j)
					{
						ci.push_back(li[j] - sbc);
						ci.push_back(li[j] + sbc);
					}
				}
				indices.push_back(ci);

				std::vector<int32_t> cl;
				std::vector<BinaryTreeNode<int32_t>*> nextLevel;
				for (int32_t j = 0; j < currentLevel.size(); ++j)
				{
					if (currentLevel[j] == nullptr)
					{
						cl.push_back(-1);
						nextLevel.push_back(nullptr);
						nextLevel.push_back(nullptr);
					}
					else
					{
						cl.push_back(currentLevel[j]->Key());
						nextLevel.push_back(currentLevel[j]->Left());
						nextLevel.push_back(currentLevel[j]->Right());
					}
						
				}
				currentLevel.clear();
				currentLevel = nextLevel;
				levels.push_back(cl);
			}

			std::vector<std::string> strings;
			// initialize
			int32_t count = std::pow(2, h) - 1;
			for (int32_t i = 0; i < count; ++i)
				strings.push_back(std::string());

			// print the keys
			for (int32_t i = 0; i < levels.size(); ++i)
			{
				std::vector<int32_t>& level = levels[i];
				std::vector<int32_t>& index = indices[i];
				int32_t t = 0;
				for (int32_t j = 0; j < level.size(); ++j)
				{
					char chs[6] = { 0 };
					itoa(level[j], chs, 10);
					std::string oStr(chs);
					std::string tStr("     ");
					for (int32_t k = oStr.size() - 1, l = 4; k >= 0; --k, --l)
						tStr[l] = oStr[k];
					if (i != 0)
						tStr[0] = '-';
					if (oStr == "-1")
					{
						tStr = "     ";
					}
					strings[index[j]].append(tStr);
					
					for (; t < index[j]; ++t)
					{
						strings[t].append("     ");
					}
					t = index[j] + 1;
				}

				// fill the left strings
				for (; t < strings.size(); ++t)
					strings[t].append("     ");

				// fill in the bars
				if (i != levels.size() - 1)
				{
					std::vector<int32_t>& nextLevel = levels[i + 1];
					for (int32_t j = 0; j < level.size(); ++j)
					{
						int32_t l = nextLevel[2 * j];
						int32_t r = nextLevel[2 * j + 1];
						int32_t dashCount = std::pow(2, h - i - 1 - 1);
						if (l != -1)
							for (int32_t k = 1; k < dashCount; ++k)
								strings[index[j] - k][5 * (i + 1) - 1] = '|';

						if (r != -1)
							for (int32_t k = 1; k < dashCount; ++k)
								strings[index[j] + k][5 * (i + 1) - 1] = '|';
					}
				}
			}

			// to file
			for (int32_t i = 0; i < strings.size(); ++i)
			{
				out << strings[i].c_str() << std::endl;
			}
		}

		void DumpWithPostorderWalk(std::ostream& out)
		{
			out << "【Postorder Walk】The tree is : \n";
			PostorderWalk(out, m_Root);
			out << std::endl << std::endl;
		}

		void DumpWithInvertorderWalk(std::ostream& out)
		{
			out << "【InvertOrder Walk】The tree is:\n";
			InvertorderWalk(out, m_Root);
			out << std::endl << std::endl;
		}

		int32_t Height()
		{
			return Height(m_Root);
		}

		int32_t Height(BinaryTreeNode<T>* node)
		{
			if (node == nullptr)
				return 0;

			int32_t left = Height(node->Left());
			int32_t right = Height(node->Right());

			return left >= right ? left + 1: right + 1;
		}

    protected:
        void InorderWalk(std::ostream& out, BinaryTreeNode<T>* node)
        {
            if (node->Left() != nullptr)
                InorderWalk(out, node->Left());
            
            out << node->Key() << ' ';

            if (node->Right() != nullptr)
                InorderWalk(out, node->Right());
        }

		void InorderWalk(std::vector<T>& list, BinaryTreeNode<T>* node)
		{
			if (node->Left() != nullptr)
				InorderWalk(list, node->Left());

			list.push_back(node->Key());

			if (node->Right() != nullptr)
				InorderWalk(list, node->Right());
		}

		void PreorderWalk(std::ostream& out, BinaryTreeNode<T>* node)
		{
			out << node->Key() << ' ';

			if (node->Left() != nullptr)
				PreorderWalk(out, node->Left());

			if (node->Right() != nullptr)
				PreorderWalk(out, node->Right());
		}

		void PostorderWalk(std::ostream& out, BinaryTreeNode<T>* node)
		{
			if (node->Left() != nullptr)
				PostorderWalk(out, node->Left());
			if (node->Right() != nullptr)
				PostorderWalk(out, node->Right());

			out << node->Key() << ' ';
		}

		void InvertorderWalk(std::ostream& out, BinaryTreeNode<T>* node)
		{
			if (node->Right() != nullptr)
				InvertorderWalk(out, node->Right());

			out << node->Key() << ' ';

			if (node->Left() != nullptr)
				InvertorderWalk(out, node->Left());
		}

		void Free(BinaryTreeNode<T>* node)
		{
			if (node->Left() != nullptr)
				Free(node->Left());

			if (node->Right() != nullptr)
				Free(node->Right());

			delete node;
		}

    protected:
        BinaryTreeNode<T>*      m_Root;
    };

    // Another binary search tree, works well but needs lots of task than the BST above
    // Red-Black tree
    /**
     * The properties of red-black tree is:
     * 1. Every node is either red or black.
     * 2. The root is black
     * 3. Every leaf is black
     * 4. If a node is red, then both its children are black.
     * 5. For each node, all simple paths from the node to descendant leaves 
     *    contain the same number of black nodes.
     * */

    // First we will define red-black node
    // which is a kind of BinarySearchNode
    template <typename T>
    class RedBlackNode
    {
    public:
        RedBlackNode() : m_Parent(nullptr), m_Left(nullptr), m_Right(nullptr), m_Key(T(0)), m_IsRed(true) {}
		RedBlackNode(T key) : m_Parent(nullptr), m_Left(nullptr), m_Right(nullptr), m_Key(key), m_IsRed(true) {}
        RedBlackNode(RedBlackNode* p, RedBlackNode* l, RedBlackNode* r, T key, bool red)
            : m_Parent(p), m_Left(l), m_Right(r), m_Key(key), m_IsRed(red) {}

    public:
        RedBlackNode* Parent()
        {
            return m_Parent;
        }
        void Parent(RedBlackNode* p)
        {
            m_Parent = p;
        }

        RedBlackNode* Left()
        {
            return m_Left;
        }
        void Left(RedBlackNode* l)
        {
            m_Left = l;
        }

        RedBlackNode* Right()
        {
            return m_Right;
        }
        void Right(RedBlackNode* r)
        {
            m_Right = r;
        }

        T Key()
        {
            return m_Key;
        }
        void Key(T key)
        {
            m_Key = key;
        }

        bool IsRed()
        {
            return m_IsRed;
        }
        void SetRed()
        {
            m_IsRed = true;
        }
        void SetBlack()
        {
            m_IsRed = false;
        }
    
    protected:
        virtual void Dump(std::ostream& out) const {}

	protected:
		RedBlackNode* m_Parent;
		RedBlackNode* m_Left;
		RedBlackNode* m_Right;
		T m_Key;
        bool m_IsRed;
    };
    template<typename T>
    class RedBlackTree
    {
    public:
        RedBlackTree()
		{
			m_Nil.SetBlack();
			m_Root = &m_Nil;
		}
        RedBlackTree(BinaryTreeNode<T>* root)
		{
			m_Nil.SetBlack();
			m_Root = &m_Nil;
		}
        virtual ~RedBlackTree() {}

    public:
        bool IsEmpty() {return m_Root == &m_Nil;}
        RedBlackNode<T>* Root() { return m_Root; }
        RedBlackNode<T>* Search(RedBlackNode<T>* p, T key)
        {
			if (p == nullptr)
				p = &m_Nil;
            if (p == &m_Nil || key == p->Key())
                return p;
            else if (key <= p->Key())
                return Search(p->Left(), key);
            else 
                return Search(p->Right(), key);
        }

        RedBlackNode<T>* Minimum(RedBlackNode<T>* p)
        {
            RedBlackNode<T>* t = p? p : &m_Nil;
            while (t->Left() != &m_Nil)
                t = t->Left();
            return t;
        }

        RedBlackNode<T>* Maximum(RedBlackNode<T>* p)
        {
            RedBlackNode<T>* t = p? p : &m_Nil;
            while (t->Right() != &m_Nil)
                t = t->Right();
            return t;
        }

        RedBlackNode<T>* Successor(RedBlackNode<T>* p)
        {
            if (p->Right() != &m_Nil)
                return Minimum(p->Right());

            RedBlackNode<T>* t1 = p;
            RedBlackNode<T>* t2 = p->Parent();
            while (t2 != &m_Nil && t2->Right() == t1)
            {
                t1 = t2;
                t2 = t2->Parent();
            }

            return t2;
        }

        RedBlackNode<T>* Predecessor(RedBlackNode<T>* p)
        {
            if (p->Left() != &m_Nil)
                return Maximum(p->Left());

            RedBlackNode<T>* t1 = p;
            RedBlackNode<T>* t2 = p->Parent();
            while (t2 != &m_Nil && t2->Left() == t1)
            {
                t1 = t2;
                t2 = t2->Parent();
            }

            return t2;
        }

        void LeftRotate(RedBlackNode<T>* p)
        {
            RedBlackNode<T>* x = p;
            // deal with p's right child's left child
            RedBlackNode<T>* y = x->Right();
            x->Right(y->Left());
            if (x->Right() != &m_Nil)
                x->Right()->Parent(x);

            // deal with p's parent
            y->Parent(x->Parent());
            if (y->Parent() == &m_Nil)
                m_Root = y;
            else if (y->Parent()->Left() == x)
                y->Parent()->Left(y);
            else 
                y->Parent()->Right(y);

            // deal with t and r
            y->Left(x);
            x->Parent(y);
        }

        void RightRotate(RedBlackNode<T>* p)
        {
            RedBlackNode<T>* y = p;
            // deal with p's left child's right child
            RedBlackNode<T>* x = y->Left();
            y->Left(x->Right());
            if (y->Left() != &m_Nil)
                y->Left()->Parent(y);

            // deal with p's parent
            x->Parent(y->Parent());
            if (x->Parent() == &m_Nil)
                m_Root = x;
            else if (x->Parent()->Left() == y)
                x->Parent()->Left(x);
            else 
                x->Parent()->Right(x);

            // deal with y and x
            x->Right(y);
            y->Parent(x);
        }

        void Transplant(RedBlackNode<T>* u, RedBlackNode<T>* v)
        {
            if (u->Parent() == &m_Nil)
                m_Root = v;
            else if (u->Parent()->Left() == u)
                u->Parent()->Left(v);
            else 
                u->Parent()->Right(v);
            v->Parent(u->Parent());
        }

        void Insert(RedBlackNode<T>* p)
        {
            RedBlackNode<T>* y = &m_Nil;
            RedBlackNode<T>* x = m_Root;
            while (x != &m_Nil)
            {
                y = x;
                if (p->Key() <= x->Key())
                    x = x->Left();
                else 
                    x = x->Right();
            }
            p->Parent(y);
            if (y == &m_Nil)
                m_Root = p;
            else if (p->Key() <= y->Key())
                y->Left(p);
            else 
                y->Right(p);
            p->Left(&m_Nil);
            p->Right(&m_Nil);
            p->SetRed();
            InsertFixup(p);
        }

        void Delete(RedBlackNode<T>* p)
        {
            RedBlackNode<T>* x = &m_Nil;
            RedBlackNode<T>* y = p;
            bool flagIsRed = y->IsRed();
            if (p->Left() == &m_Nil) // if p doesn't have a left child
            {
                x = p->Right(); // Now, x will replcace the position of p
                Transplant(p, p->Right());
            }
            else if (p->Right() == &m_Nil)
            {
                x = p->Left(); // Also, x will replcace the position of p
                Transplant(p, p->Left());
            }
            else 
            {
                y = Minimum(p->Right()); // Now, y is successor of p
                flagIsRed = y->IsRed(); // Set y's color to be the flag. Cause we will set y to be p's color later.
                x = y->Right(); // y can have right child but it can't have left child
                if (y->Parent() == p) // This may happen when y is p's right child
                    x->Parent(y); // y will replace p, so x is still y's child
                else 
                { // Or, x will replace y
                    Transplant(y, x);
                    y->Right(p->Right()); // y obtains p's right
                    y->Right()->Parent(y);
                }

                Transplant(p, y); // replace p with y
                y->Left(p->Left()); // y obtains p's left
                y->Left()->Parent(y);
                if (p->IsRed()) y->SetRed(); // set y's color to p's color
                else y->SetBlack();
            }

            if (!flagIsRed) // There are only two cases which will cause properites be violated
                            // First, p has only one child, and it's a black node. We deleted it so we are in trouble.
                            // Second, p has two children, but we move a black node (which is y) to p's position and
                            //         put x into the orgin place of y.
                            // Conclude, if we delete or move a black node, we need to maintain the propreties.
                DeleteFixup(x);
        }

		void DumpToInorderList(std::vector<T>& list)
		{
			InorderWalk(list, m_Root);
		}

		void DumpWithInorderWalk(std::ostream& out)
		{
			out << "【Inorder Walk】 The tree is: \n";
			InorderWalk(out, m_Root);
			out << std::endl << std::endl;
		}

		void DumpWithPreorderWalk(std::ostream& out)
		{
			out << "【Preorder Walk】The tree is : \n";
			PreorderWalk(out, m_Root);
			out << std::endl << std::endl;
		}

		void DumpStructureToFile(std::fstream& out)
		{
			int32_t h = Height();
			assert(h != 0);
			if (h == 0)
				return;

			std::vector<std::vector<int32_t>> levels;
			std::vector<std::vector<int32_t>> indices;
			std::vector<RedBlackNode<int32_t>*> currentLevel;
			currentLevel.push_back(m_Root);
			for (int32_t i = 1; i <= h; ++i)
			{
				// The indices of current level
				std::vector<int32_t> ci;
				std::vector<int32_t> li;
				if (indices.size() > 0)
					li = indices[indices.size() - 1];
				if (li.size() == 0) // We are now dealing with root level.
					ci.push_back(std::pow(2, h - 1) - 1);
				else
				{
					int32_t sbc = std::pow(2, h - i);// side bar count of last level 
					for (int32_t j = 0; j < li.size(); ++j)
					{
						ci.push_back(li[j] - sbc);
						ci.push_back(li[j] + sbc);
					}
				}
				indices.push_back(ci);

				std::vector<int32_t> cl;
				std::vector<RedBlackNode<int32_t>*> nextLevel;
				for (int32_t j = 0; j < currentLevel.size(); ++j)
				{
					if (currentLevel[j] == &m_Nil)
					{
						cl.push_back(-1);
						nextLevel.push_back(&m_Nil);
						nextLevel.push_back(&m_Nil);
					}
					else
					{
						cl.push_back(currentLevel[j]->Key());
						nextLevel.push_back(currentLevel[j]->Left());
						nextLevel.push_back(currentLevel[j]->Right());
					}

				}
				currentLevel.clear();
				currentLevel = nextLevel;
				levels.push_back(cl);
			}

			std::vector<std::string> strings;
			// initialize
			int32_t count = std::pow(2, h) - 1;
			for (int32_t i = 0; i < count; ++i)
				strings.push_back(std::string());

			// print the keys
			for (int32_t i = 0; i < levels.size(); ++i)
			{
				std::vector<int32_t>& level = levels[i];
				std::vector<int32_t>& index = indices[i];
				int32_t t = 0;
				for (int32_t j = 0; j < level.size(); ++j)
				{
					char chs[6] = { 0 };
					itoa(level[j], chs, 10);
					std::string oStr(chs);
					std::string tStr("     ");
					for (int32_t k = oStr.size() - 1, l = 4; k >= 0; --k, --l)
						tStr[l] = oStr[k];
					if (i != 0)
						tStr[0] = '-';
					if (oStr == "-1")
					{
						tStr = "     ";
					}
					strings[index[j]].append(tStr);

					for (; t < index[j]; ++t)
					{
						strings[t].append("     ");
					}
					t = index[j] + 1;
				}

				// fill the left strings
				for (; t < strings.size(); ++t)
					strings[t].append("     ");

				// fill in the bars
				if (i != levels.size() - 1)
				{
					std::vector<int32_t>& nextLevel = levels[i + 1];
					for (int32_t j = 0; j < level.size(); ++j)
					{
						int32_t l = nextLevel[2 * j];
						int32_t r = nextLevel[2 * j + 1];
						int32_t dashCount = std::pow(2, h - i - 1 - 1);
						if (l != -1)
							for (int32_t k = 1; k < dashCount; ++k)
								strings[index[j] - k][5 * (i + 1) - 1] = '|';

						if (r != -1)
							for (int32_t k = 1; k < dashCount; ++k)
								strings[index[j] + k][5 * (i + 1) - 1] = '|';
					}
				}
			}

			// to file
			for (int32_t i = 0; i < strings.size(); ++i)
			{
				out << strings[i].c_str() << std::endl;
			}
		}

        void DumpRedBlackLevelsToFile(std::fstream& out)
        {
			int32_t h = Height();
			assert(h != 0);
			if (h == 0)
				return;

			std::vector<std::vector<int32_t>> levels;
			std::vector<std::vector<int32_t>> indices;
			std::vector<RedBlackNode<int32_t>*> currentLevel;
			currentLevel.push_back(m_Root);
			for (int32_t i = 1; i <= h; ++i)
			{
				// The indices of current level
				std::vector<int32_t> ci;
				std::vector<int32_t> li;
				if (indices.size() > 0)
					li = indices[indices.size() - 1];
				if (li.size() == 0) // We are now dealing with root level.
					ci.push_back(std::pow(2, h - 1) - 1);
				else
				{
					int32_t sbc = std::pow(2, h - i);// side bar count of last level 
					for (int32_t j = 0; j < li.size(); ++j)
					{
						ci.push_back(li[j] - sbc);
						ci.push_back(li[j] + sbc);
					}
				}
				indices.push_back(ci);

				std::vector<int32_t> cl;
				std::vector<RedBlackNode<int32_t>*> nextLevel;
				for (int32_t j = 0; j < currentLevel.size(); ++j)
				{
					if (currentLevel[j] == &m_Nil)
					{
						cl.push_back(-1);
						nextLevel.push_back(&m_Nil);
						nextLevel.push_back(&m_Nil);
					}
					else
					{
						cl.push_back(currentLevel[j]->IsRed()? 1 : 0);
						nextLevel.push_back(currentLevel[j]->Left());
						nextLevel.push_back(currentLevel[j]->Right());
					}

				}
				currentLevel.clear();
				currentLevel = nextLevel;
				levels.push_back(cl);
			}

			std::vector<std::string> strings;
			// initialize
			int32_t count = std::pow(2, h) - 1;
			for (int32_t i = 0; i < count; ++i)
				strings.push_back(std::string());

			// print the keys
			for (int32_t i = 0; i < levels.size(); ++i)
			{
				std::vector<int32_t>& level = levels[i];
				std::vector<int32_t>& index = indices[i];
				int32_t t = 0;
				for (int32_t j = 0; j < level.size(); ++j)
				{
					std::string tStr("  ");
					if (i != 0 && level[j] != -1)
						tStr[0] = '-';
                    if (level[j] == 0)
                        tStr[1] = 'b';
                    else if (level[j] == 1)
                        tStr[1] = 'r';
                    else 
                        tStr[1] = ' ';
					strings[index[j]].append(tStr);

					for (; t < index[j]; ++t)
					{
						strings[t].append("  ");
					}
					t = index[j] + 1;
				}

				// fill the left strings
				for (; t < strings.size(); ++t)
					strings[t].append("  ");

				// fill in the bars
				if (i != levels.size() - 1)
				{
					std::vector<int32_t>& nextLevel = levels[i + 1];
					for (int32_t j = 0; j < level.size(); ++j)
					{
						int32_t l = nextLevel[2 * j];
						int32_t r = nextLevel[2 * j + 1];
						int32_t dashCount = std::pow(2, h - i - 1 - 1);
						if (l != -1)
							for (int32_t k = 1; k < dashCount; ++k)
								strings[index[j] - k][2 * (i + 1) - 1] = '|';

						if (r != -1)
							for (int32_t k = 1; k < dashCount; ++k)
								strings[index[j] + k][2 * (i + 1) - 1] = '|';
					}
				}
			}

			// to file
			for (int32_t i = 0; i < strings.size(); ++i)
			{
				out << strings[i].c_str() << std::endl;
			}
        }

		void DumpWithPostorderWalk(std::ostream& out)
		{
			out << "【Postorder Walk】The tree is : \n";
			PostorderWalk(out, m_Root);
			out << std::endl << std::endl;
		}

		void DumpWithInvertorderWalk(std::ostream& out)
		{
			out << "【InvertOrder Walk】The tree is:\n";
			InvertorderWalk(out, m_Root);
			out << std::endl << std::endl;
		}

		int32_t Height()
		{
			return Height(m_Root);
		}

		int32_t Height(RedBlackNode<T>* node)
		{
			if (node == nullptr)
				return 0;

			int32_t left = Height(node->Left());
			int32_t right = Height(node->Right());

			return left >= right ? left + 1 : right + 1;
		}

    protected:
        void InsertFixup(RedBlackNode<T>* p)
        {
            RedBlackNode<T>* t = p;
            while (t->Parent()->IsRed())
            {
                if (t->Parent() == t->Parent()->Parent()->Left())
                {
                    RedBlackNode<T>* y = t->Parent()->Parent()->Right();
                    if (y->IsRed()) // case 1: if my 'uncle' node is red
                    {
                        t->Parent()->SetBlack();
                        y->SetBlack();
                        t->Parent()->Parent()->SetRed();
                        t = t->Parent()->Parent();
                    }
                    else if (t == t->Parent()->Right()) // case 2: my 'uncle' is black, I am the right child of my parent
                    {
                        t = t->Parent();
                        LeftRotate(t);
                    }
                    else // case 3: my 'uncle' is black, I am the left child of my parent
                    {
                        t->Parent()->SetBlack();
                        t->Parent()->Parent()->SetRed();
                        RightRotate(t->Parent()->Parent());
                    }
                }
                else 
                {
                    RedBlackNode<T>* y = t->Parent()->Parent()->Left();
                    if (y->IsRed()) // case 1: if our 'uncle' node is red
                    {
                        t->Parent()->SetBlack();
                        y->SetBlack();
                        t->Parent()->Parent()->SetRed();
                        t = t->Parent()->Parent();
                    }
                    else if (t == t->Parent()->Left()) // case 2: my 'uncle' is black, I am the left child
                    {
                        t = t->Parent();
                        RightRotate(t);
                    }
                    else // case 3: my 'uncle' is black, I am the right child
                    {
                        t->Parent()->SetBlack();
                        t->Parent()->Parent()->SetRed();
                        LeftRotate(t->Parent()->Parent());
                    }
                }
            }

			m_Root->SetBlack();
        }

        void DeleteFixup(RedBlackNode<T>* x)
        {
            while (x != m_Root && !x->IsRed())
            {
                if (x == x->Parent()->Left())
                {
                    RedBlackNode<T>* w = x->Parent()->Right();
                    if (w->IsRed())
                    {   // case 1
                        w->SetBlack();
                        x->Parent()->SetRed();
                        LeftRotate(x->Parent());
                        w = x->Parent()->Right();
                    }

                    if (!w->Left()->IsRed() && !w->Right()->IsRed())
                    {   // case 2
                        w->SetRed();
                        x = x->Parent();
                    }
                    else if (!w->Right()->IsRed())
                    {   // case 3
                        w->Left()->SetBlack();
                        w->SetRed();
                        RightRotate(w);
                        w = x->Parent()->Right();
                    }
                    else 
                    {   // case 4
                        if (x->Parent()->IsRed()) w->SetRed();
                        else w->SetBlack();
                        x->Parent()->SetBlack();
                        w->Right()->SetBlack();
                        LeftRotate(x->Parent());
                        x = m_Root;
                    }
                }
                else 
                {
                    RedBlackNode<T>* w = x->Parent()->Left();
                    if (w->IsRed())
                    {   // case 1
                        w->SetBlack();
                        x->Parent()->SetRed();
                        RightRotate(x->Parent());
                        w = x->Parent()->Left();
                    }

                    if (!w->Right()->IsRed() && !w->Left()->IsRed())
                    {   // case 2
                        w->SetRed();
                        x = x->Parent();
                    }
                    else if (!w->Left()->IsRed())
                    {   // case 3
                        w->Right()->SetBlack();
                        w->SetRed();
                        LeftRotate(w);
                        w = x->Parent()->Left();
                    }
                    else 
                    {   // case 4
                        if (x->Parent()->IsRed()) w->SetRed();
                        else w->SetBlack();
                        x->Parent()->SetBlack();
                        w->Right()->SetBlack();
                        RightRotate(x->Parent());
                        x = m_Root;
                    }
                }
            }

            x->SetBlack();
        }

		void InorderWalk(std::ostream& out, RedBlackNode<T>* node)
		{
			if (node->Left() != &m_Nil)
				InorderWalk(out, node->Left());

			out << node->Key() << ' ';

			if (node->Right() != &m_Nil)
				InorderWalk(out, node->Right());
		}

		void InorderWalk(std::vector<T>& list, RedBlackNode<T>* node)
		{
			if (node->Left() != &m_Nil)
				InorderWalk(list, node->Left());

			list.push_back(node->Key());

			if (node->Right() != &m_Nil)
				InorderWalk(list, node->Right());
		}

		void PreorderWalk(std::ostream& out, RedBlackNode<T>* node)
		{
			out << node->Key() << ' ';

			if (node->Left() != &m_Nil)
				PreorderWalk(out, node->Left());

			if (node->Right() != &m_Nil)
				PreorderWalk(out, node->Right());
		}

		void PostorderWalk(std::ostream& out, RedBlackNode<T>* node)
		{
			if (node->Left() != &m_Nil)
				PostorderWalk(out, node->Left());
			if (node->Right() != &m_Nil)
				PostorderWalk(out, node->Right());

			out << node->Key() << ' ';
		}

		void InvertorderWalk(std::ostream& out, RedBlackNode<T>* node)
		{
			if (node->Right() != &m_Nil)
				InvertorderWalk(out, node->Right());

			out << node->Key() << ' ';

			if (node->Left() != &m_Nil)
				InvertorderWalk(out, node->Left());
		}

		void Free(RedBlackNode<T>* node)
		{
			if (node->Left() != &m_Nil)
				Free(node->Left());

			if (node->Right() != &m_Nil)
				Free(node->Right());

			delete node;
		}

    protected:
        RedBlackNode<T>*	m_Root;
		RedBlackNode<T>		m_Nil;
    };
}