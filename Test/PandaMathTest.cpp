#include <iostream>
#include <assert.h>
#include <iomanip>
#include "Math/PandaMath.hpp"
#include <fstream>
#include <sstream>
#include <iterator>
#include "HighPrecisionTimer.hpp"
#include "Math/Tree.hpp"
#include "Math/SortAlgos.hpp"
#include <ctime>	// std::time
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace Panda;
using namespace std;

namespace Panda
{
	Handness g_ViewHandness = Handness::kHandnessRight;
	DepthClipSpace g_DepthClipSpace = DepthClipSpace::kDepthClipNegativeOneToOne;
}

template<typename T>
ostream& Output (ostream& out, int32_t cut, int32_t count, T in[])
{
    for (int i = 0; i < count; ++i)
    {
		if (i != 0 && (i % cut == 0))
			out << endl;
        out << fixed << setprecision(6) <<  in[i] << ",";
    }

	out << endl;
    return out;
}

void VectorTest()
{
	Vector2Df x({ 55.3f, 22.1f });
    cout << "Vector2Df: " << x << endl << endl;

    Vector3Df a({1.0f, 2.0f, 3.0f});
    Vector3Df b({5.0f, 6.0f, 7.0f});

    cout << "Vector a : " << a << endl;
    cout << "Vector b : " << b << endl;

    Vector3Df c;
    c = CrossProduct(a, b);
    cout<< "Cross Product of vec a and b is: " << c << endl;
    
    float d;
    d = DotProduct(a, b);
    cout << "Dot Product of vec a and b is: " << d << endl << endl;

    Vector4Df e({-3.0f, 3.0f, 6.0f, 1.0f});
    Vector4Df f({2.0f, 0.0f, -0.7f, 0.0f});
    cout << "Vector e: " << e << endl;
    cout << "Vector f: " << f << endl;

    Vector4Df g = e + f;
    cout << "vec e + vec f: " << g << endl;
    g = e - f;
    cout << "vec e - vec f: "<< g << endl;

    g = Normalize(g);
    cout << "g normalized is " << g << endl << endl;
}

void MatrixTest()
{
	Matrix4f m1;
	m1.SetIdentity();

	cout << "Identity Matrix : " << m1 << endl;

	Matrix4f mEu;
	float yaw = 0.2f, pitch = 0.3f, roll = 0.4f;
	MatrixRotationYawPitchRoll(mEu, yaw, pitch, roll);
	cout << "Matrix of yaw(" << yaw << ") pitch(" << pitch << ") roll(" << roll << "):" << mEu << endl;

	Matrix4f ry;
	float angle = PI / 2.0f;
	MatrixRotationY(ry, angle);
	cout << "Matrix of Rotation on Y(angle = " << angle << "): " << ry << endl;

	Matrix4f rz;
	MatrixRotationZ(rz, angle);
	cout << "Matrix of Rotation on Z(angle = " << angle << "): " << rz << endl;

	float x = 5.0f, y = 6.5f, z = -7.0f;
	Matrix4f translate;
	MatrixTranslation(translate, x, y, z);
	cout << "Matrix of Translation on X(" << x << ")Y" << y << ")Z(" << z << "):";
	cout << translate << endl;

	cout << "Matrix multiplication: Rotation Y * Rotation Z * Translation on X: ";
	Matrix4f transform = m1 * ry * rz * translate;
	cout << transform << endl;

	Vector3Df v({ 1.0f, 0.0f, 0.0f });

	Vector3Df v1 = v;
	cout << "Vector : " << v1;
	cout << "Transfrom by Rotation Y Matrix: ";
	cout << ry;
	TransformCoord(v1, ry);
	cout << "Now the vector becomes: " << v1;
	cout << endl;

	v1 = v;
	cout << "Vector: " << v1;
	cout << "Tranform by Rotation Z Matrix: ";
	cout << rz;
	TransformCoord(v1, rz);
	cout << "Now the vector becomes: " << v1;
	cout << endl;

	v1 = v;
	cout << "Vector: " << v1;
	cout << "Transform by Translation Matrix: ";
	cout << translate;
	TransformCoord(v1, translate);
	cout << "Now the vector becomes: " << v1;
	cout << endl;

	Vector3Df position({ 0.0f, 0.0f, -5.0f }), lookAt({ 0.0f, 0.0f, 0.0f }), up({ 0.0f, 1.0f, 0.0f });
	Matrix4f view;
	BuildViewMatrix(view, position, lookAt, up);
	cout << "View Matrix with positoin(" << position << ") lookAt(" << lookAt << ") up(" << up << "):";
	cout << view << endl;

	float fov = PI / 2.0f, aspect = 16.0f / 9.0f, near = 1.0f, far = 100.0f;
	Matrix4f perspective;
	BuildPerspectiveFovLHMatrix(perspective, fov, aspect, near, far);
	cout << "(Left-Handed) Perspective Matrix with fov(" << fov << ") aspect(" << aspect << ") near ... far(" << near << "..." << far << "):";
	cout << perspective << endl;
	BuildPerspectiveFovRHMatrix(perspective, fov, aspect, near, far);
	cout << "(Right-Handed) Perspective Matrix with fov(" << fov << ") aspect(" << aspect << ") near ... far(" << near << "..." << far << "):";
	cout << perspective << endl;

	Matrix4f mvp = perspective * view;
	cout << "MVP: " << mvp << endl;

	Matrix3f mat3(
		{ 2.0f, 1.0f, 10.0f,
		1.0f, 2.0f, 3.0f,
		2.0f, 3.0f, 4.0f }
	);
	InverseMatrix(mat3, mat3);
	cout << "Inverse mat3 : " << mat3 << endl;

    Matrix4f invertable(
		{ 1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		13.0f, 14.0f, 15.0f, 1.0f });
    cout << "Known Invertable Matrix: " << invertable;
    assert(InverseMatrix(invertable, invertable));
    cout << "Inverse of Matrix: " << invertable;

    Matrix4f nonInvertable({
         1.0f,  2.0f,  3.0f,  4.0f,
         5.0f,  6.0f,  7.0f,  8.0f,
         9.0f, 10.0f, 11.0f, 12.0f,
        13.0f, 14.0f, 15.0f, 16.0f
		});
    cout << "Know sigular Matrix: " << nonInvertable;
    assert(!InverseMatrix(nonInvertable, nonInvertable));
    cout << "InvertMatrix4f returns false. " << endl;
}

void TestDCTs()
{
    float pixelBlock[64] = {
        -76, -73, -67, -62, -58, -67, -64, -55,
        -65, -69, -73, -38, -19, -43, -59, -56,
        -66, -69, -60, -15,  16, -24, -62, -55,
        -65, -70, -57, - 6,  26, -22, -58, -59,
        -61, -67, -60, -24, - 2, -40, -60, -58,
        -49, -63, -68, -58, -51, -60, -70, -53,
        -43, -57, -64, -69, -73, -67, -63, -45,
        -41, -49, -59, -60, -63, -52, -50, -34
    };

    cout << "A 8x8 int pixel block: " << endl;
    Output(cout, 8, 64, pixelBlock);
    float dct[64] = {0.0};
    DCT8x8(pixelBlock, dct);
	cout << endl;
    cout << "After DCTII : " << endl;
    Output(cout, 8, 64, dct);
	cout << endl;
	cout << endl;

	float pixelBlock2[64] = {
		-416, -33, -60,  32,  48, -40,   0,   0,
		   0, -24, -56,  19,  26,   0,   0,   0,
		 -42,  13,  80, -24, -40,   0,   0,   0,
		 -42,  17,  44, -29,   0,   0,   0,   0,
		  18,   0,   0,   0,   0,   0,   0,   0,
		   0,   0,   0,   0,   0,   0,   0,   0,
		   0,   0,   0,   0,   0,   0,   0,   0,
		   0,   0,   0,   0,   0,   0,   0,   0
	};

	cout << "A 8x8 int pixel block 2: " << endl;
	Output(cout, 8, 64, pixelBlock2);
	float dct2[64] = { 0.0f };
	IDCT8x8(pixelBlock2, dct2);
	cout << endl;
	cout << "After DCTIII : " << endl;
	Output(cout, 8, 64, dct2);
	cout << endl;
}

void TestQuaternion()
{
	float root2 = sqrtf(2.0f);
	float halfRoot2 = root2 / 2.0f;
	Vector4Df p({ 1.0f, 0.0f, 0.0f, 1.0f });	// origin point <1, 0, 0>
	Vector4Df q({ 0.0f, 0.0f, halfRoot2, halfRoot2 });	 // transform along z-axis rotate 90 degrees
	Vector4Df q_({ 0.0f, 0.0f, -halfRoot2, halfRoot2 }); // Conjugate of q
	Vector4Df p1 = MulByElement(MulByElement(q, p), q_);

	cout << "p1 = " << p1 << endl;

	Vector4Df r({ halfRoot2, 0, 0, halfRoot2 }); // along x-axis rotate 90 degrees
	Vector4Df r_({ -halfRoot2, 0, 0, halfRoot2 });	 // conjugate of r
	Vector4Df p2 = MulByElement(MulByElement(r, p1), r_);
	cout << "p2 = " << p2 << endl;

	Vector4Df w = MulByElement(r, q);
	Vector4Df w_ = MulByElement(q_, r_);
	Vector4Df p3 = MulByElement(MulByElement(w, p), w_);
	cout << "w = " << w << endl;
	cout << "w_ = " << w_ << endl;
	cout << "w * p = " << MulByElement(w, p) << endl;
	cout << "p3 = " << p3 << endl;

	Vector4Df v({ 0.0f, -halfRoot2, 0.0f, halfRoot2 });
	Vector4Df v_({ 0.0f, halfRoot2, 0.0f, halfRoot2 });
	Vector4Df p4 = MulByElement(MulByElement(v, p), v_);
	cout << "p4 = " << p4 << endl;
}

void TestHeap()
{
	std::vector<int32_t> arr = {20, 39, 2, 93, 9, 2, 9, 63, 9, 24, 11, 19, 39, 48, 46, 92, 24};
	std::vector<int32_t> sortedArray = HeapSortDecreasing(arr);
	for (int32_t i = 0; i < sortedArray.size(); ++i)
		cout << sortedArray[i] << " ";
}

void TestQuickSort()
{
	srand(time(0));
	std::vector<int32_t> arr = { 20, 39, 2, 93, 9, 2, 9, 63, 9, 24, 11, 19, 39, 48, 46, 92, 24 };
	cout << "Before quick sorting : " << endl;
	for (int32_t i = 0; i < arr.size(); ++i)
		cout << arr[i] << " ";
	cout << endl << endl << endl;
	QuickSortDecreasing(arr, 0, arr.size() - 1);
	cout << "After decreasing quick sorting: " << endl;
	for (int32_t i = 0; i < arr.size(); ++i)
		cout << arr[i] << " ";
	cout << endl << endl << endl;
	
	QuickSortIncreasing(arr, 0, arr.size() - 1);
	cout << "After increasing quick sorting: " << endl;
	for (int32_t i = 0; i < arr.size(); ++i)
		cout << arr[i] << " ";
	cout << endl;
}

void Generate1MillionNumbers()
{
	srand(time(0));
	fstream file;
	file.open("numbers.txt", std::ios_base::out);
	for (int32_t i = 0; i < 1000000; ++i)
		file << rand() % 100000 << ' ';
	file.close();
}

void TestQuickSort1()
{
	HighPrecisionTimer timer;

	cout << "Reading the file" << endl;
	ifstream inFile;
	inFile.open("numbers.txt"); // default open mode: as chars
	string str;
	getline(inFile, str);
	inFile.close();
	
	// split the words 
	std::vector<std::string> vec;
	istringstream iss(str);
	copy(istream_iterator<string>(iss),
		istream_iterator<string>(),
		back_inserter(vec));
	
	// into integers
	std::vector<int32_t> ints;
	for (string& str : vec)
	{
		ints.push_back(atoi(str.c_str()));
	}

	// Start Sorting
	cout << "Start Sorting ..." << endl;
	timer.Start();
	QuickSortIncreasing(ints, 0, ints.size() - 1);
	timer.Stop();
	cout << "End Sorting ..." << endl;
	cout << "Cost time " << timer.TotalTime() / 1000.0f << "seconds" << endl;

	cout << endl << endl;
	cout << "Start writing back to file" << endl;
	ofstream oFile("numbers_.txt");
	for (int32_t i = 0; i < ints.size(); ++i)
		oFile << ints[i] << ' ';
	oFile.close();
}

void TestCountingSort()
{
	HighPrecisionTimer timer;

	cout << "Reading the file" << endl;
	ifstream inFile;
	inFile.open("numbers.txt"); // default open mode: as chars
	string str;
	getline(inFile, str);
	inFile.close();

	// split the words 
	std::vector<std::string> vec;
	istringstream iss(str);
	copy(istream_iterator<string>(iss),
		istream_iterator<string>(),
		back_inserter(vec));

	// into integers
	std::vector<int32_t> ints;
	for (string& str : vec)
	{
		ints.push_back(atoi(str.c_str()));
	}

	// Start Sorting
	cout << "Start Sorting ..." << endl;
	timer.Start();
	std::vector<int32_t> out;
	CountingSortIncreasing(ints, 32767, out);
	timer.Stop();
	cout << "End Sorting ..." << endl;
	cout << "Cost time " << timer.TotalTime() / 1000.0f << "seconds" << endl;

	cout << endl << endl;
	cout << "Start writing back to file" << endl;
	ofstream oFile("numbers_C.txt");
	for (int32_t i = 0; i < ints.size(); ++i)
		oFile << out[i] << ' ';
	oFile.close();
}

void UseFindingSecondSmallestValue()
{
	std::vector<int32_t> arr = { 20, 39, 2, 93, 9, 2, 9, 63, 9, 24, 11, 19, 39, 48, 46, 92, 24 };
	for (int32_t i = 0; i < arr.size(); ++i)
	{
		if (i == 6)
			int32_t x = 1;
		int32_t secondSmallest = RandomizeSelect(arr, 0, arr.size() - 1, i);
		QuickSortIncreasing(arr, 0, arr.size() - 1);
		assert(secondSmallest == arr[i]);
	}
}

void TestBST()
{
	srand(time(0));
	
	for (int32_t count = 0; count < 5000; ++count)
	{
		BinarySearchTree<int32_t> bst;
		for (int32_t i = 0; i < 5000; ++i)
		{
			BinaryTreeNode<int32_t>* node = new BinaryTreeNode<int32_t>(rand() % 10000);
			bst.Insert(node);
		}
		//bst.DumpWithInvertorderWalk(cout);
		//cout << "Height = " << bst.Height() << endl;

		for (int32_t i = 0; i < 100; ++i)
		{
			std::vector<int32_t> list;
			bst.DumpToInorderList(list);
			for (int32_t j = 0; j < 4999 - i; ++j)
				assert(list[j] <= list[j + 1]);
			cout << "Assert succeed! count = " << count << ", i = " << i << endl << endl;

			bst.Delete(bst.Root());
		}

		bst.Free();
	}

	//bst.DumpWithPreorderWalk(cout);
	//bst.DumpWithPostorderWalk(cout);
}

void WatchStructureOfBST()
{
	srand(time(0));

	BinarySearchTree<int32_t> bst;
	for (int32_t i = 0; i < 25; ++i)
	{
		BinaryTreeNode<int32_t>* node = new BinaryTreeNode<int32_t>(rand() % 10000);
		bst.Insert(node);
	}

	fstream oFile("1.txt");
	oFile.open("1.txt", std::ios::out);
	bst.DumpStructureToFile(oFile);
	oFile.close();
}

void TestRedBlackTree()
{
	srand(time(0));

	RedBlackTree<int32_t> rbt;
	for (int32_t i = 0; i < 500; ++i)
	{
		RedBlackNode<int32_t>* node = new RedBlackNode<int32_t>(rand() % 10000);
		node->SetRed();
		rbt.Insert(node);
	}

	for (int32_t i = 0; i < 100; ++i)
	{
		std::vector<int32_t> list;
		rbt.DumpToInorderList(list);
		for (int32_t j = 0; j < 499 - i; ++j)
			assert(list[j] <= list[j + 1]);
		cout << "Assert succeed! i = " << i << endl << endl;

		rbt.Delete(rbt.Root());
	}

	//rbt.DumpWithInorderWalk(cout);
	cout << "Height = " << rbt.Height() << endl;
}

void WatchStructureOfRBT()
{
	srand(time(0));

	RedBlackTree<int32_t> rbt;
	//for (int32_t i = 0; i < 25; ++i)
	//{
	//	RedBlackNode<int32_t>* node = new RedBlackNode<int32_t>(rand() % 10000);
	//	rbt.Insert(node);
	//}
	std::vector<int32_t> numbers({ 41, 38, 31, 12, 19, 8 });
	for (int32_t i = 0; i < numbers.size(); ++i)
	{
		RedBlackNode<int32_t>* node = new RedBlackNode<int32_t>(numbers[i]);
		rbt.Insert(node);
	}

	fstream oFile("2.txt");
	oFile.open("2.txt", std::ios::out);
	rbt.DumpStructureToFile(oFile);
	oFile.close();
}

void WatchRBStructureOfRBT()
{
	srand(time(0));

	RedBlackTree<int32_t> rbt;
	//for (int32_t i = 0; i < 25; ++i)
	//{
	//	RedBlackNode<int32_t>* node = new RedBlackNode<int32_t>(rand() % 10000);
	//	rbt.Insert(node);
	//}
	std::vector<int32_t> numbers({ 41, 38, 31, 12, 19, 8 });
	for (int32_t i = 0; i < numbers.size(); ++i)
	{
		RedBlackNode<int32_t>* node = new RedBlackNode<int32_t>(numbers[i]);
		rbt.Insert(node);
	}

	fstream oFile("3.txt");
	oFile.open("3.txt", std::ios::out);
	rbt.DumpRedBlackLevelsToFile(oFile);
	oFile.close();
}

void TestInsertSorting()
{
	Generate1MillionNumbers();

	HighPrecisionTimer timer;

	cout << "Reading the file" << endl;
	ifstream inFile;
	inFile.open("numbers.txt"); // default open mode: as chars
	string str;
	getline(inFile, str);
	inFile.close();

	// split the words 
	std::vector<std::string> vec;
	istringstream iss(str);
	copy(istream_iterator<string>(iss),
		istream_iterator<string>(),
		back_inserter(vec));

	// into integers
	std::vector<int32_t> ints;
	for (string& str : vec)
	{
		ints.push_back(atoi(str.c_str()));
	}

	// Start Sorting
	cout << "Start Sorting ..." << endl;
	timer.Start();
	InsertSortingDecreasing(ints);
	timer.Stop();
	cout << "End Sorting ..." << endl;
	cout << "Cost time " << timer.TotalTime() / 1000.0f << "seconds" << endl;

	cout << endl << endl;
	cout << "Start writing back to file" << endl;
	ofstream oFile("numbers_.txt");
	for (int32_t i = 0; i < ints.size(); ++i)
		oFile << ints[i] << ' ';
	oFile.close();
}

void TestAssignmentTime()
{
	int32_t a = 12;
	int32_t b = 17;
	int32_t* p1 = nullptr;
	int32_t* p2 = nullptr;
	int32_t* p3 = nullptr;
	p1 = &a;

	HighPrecisionTimer timer;
	timer.Start();
	for (int32_t i = 0; i < 500000; ++i)
	{
		p2 = p1;
		p3 = p1;
	}
	timer.Stop();
	cout << "The total time is : " << timer.TotalTime() << "ms" << endl;

	std::shared_ptr<int32_t> sp(new int32_t(12));
	std::shared_ptr<int32_t> sp1;
	std::shared_ptr<int32_t> sp2;
	timer.Start();
	for (int32_t i = 0; i < 500000; ++i)
	{
		sp1 = sp;
		sp2 = sp;
	}
	timer.Stop();
	cout << "The total time for shared_ptr is : " << timer.TotalTime() << "ms" << endl;
}

// This is a typedef for a random number generator
// Try boost::mt19937 or boost::ecuyer1998 instead of boost::minstd_rand
//typedef boost::minstd_rand base_generator_type;
typedef boost::random::mt19937 base_generator_type;

// This is a reproducible simulation experiment. 
void experiment(base_generator_type& generator)
{
	// Define a uniform random number distribution of integer values between
	// 1 and 6 inclusive
	typedef boost::uniform_int<> distribution_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	gen_type die_gen(generator, distribution_type(1, 6));

	// If you want to see an STL iterator interface, use iterator_adaptors.hpp
	boost::generator_iterator<gen_type> die(&die_gen);
	for (int32_t i = 0; i < 10; ++i)
		std::cout << *die++ << " ";
	std::cout << "\n";
}

// Test random
void TestRandom()
{
	// Define a random number generator and initialize it with a reproducible seed.
	base_generator_type generator(42);

	std::cout << "10 samples of a uniform distrubtuion in [0..1):\n";

	// Define a uniform random number distribution which produces "double"
	// values between 0 and 1 (0 inclusive, 1 exclusive).
	boost::uniform_real<> uni_dist(0, 1);
	boost::variate_generator<base_generator_type&, boost::uniform_real<>> uni(generator, uni_dist);

	std::cout.setf(std::ios::fixed);
	// You can now retrieve random numbers from that distribution by means
	// of a STL Generator interface, i.e. calling a generator as a zero-
	// argument function.
	for (int32_t i = 0; i < 10; ++i)
		std::cout << uni() << '\n';

	/*
	 * Change seed to something else.
	 * 
	 * (Caveat: a warning that particular things need to be considered before something can be done)
	 * (ensue: to happen after or as a result of another event)
	 * Caveat: std::time(0) is not a very good truly-random seed. When 
	 * called in rapid succession, it could return the same values and 
	 * thus the same random number sequences could ensue. If not the same
	 * values are returned, the values differ only slightly in the lowest bits.
	 * A linear congruential generator with a small factor
	 * wrapped in a uniform_smallint(see experiment) will produce the same 
	 * values for the first few iterations. This is because uniform _smallint
	 * takes only the highest bits of the generator, and the generator itself
	 * needs a few iterations to spread the initial enctropy form the lowest bits
	 * to the whole state
	*/
	generator.seed(static_cast<unsigned int>(std::time(0)));
	
	std::cout << "\nexperiment: roll a die 10 times:\n";

	// You can save a generator's state by copy construction.
	base_generator_type saved_generator = generator;

	// When calling other funcitons which take a generator or distribution
	// as a parameter, make sure to always call by reference (or pointer).
	// Calling by value invokes the copy constructor, which means that the 
	// sequence of random numbers at the caller is disconnected from the
	// sequence at the callee.
	experiment(generator);

	std::cout << "redo the experiment to verify it:\n";
	experiment(saved_generator);

	// After that, both generators are equivalent
	assert(generator == saved_generator);

	// as a degenerate case, you can set min = max for uniform_int
	boost::uniform_int<> degen_dist(4, 4);
	boost::variate_generator<base_generator_type&, boost::uniform_int<>> deg(generator, degen_dist);
	std::cout << deg() << " " << deg() << " " << deg() << std::endl;

	{
		// YOu can save the generator state for future use. You can read the state
		// back in at any later time using operator>>
		std::ofstream file("rng.saved", std::ofstream::trunc);
		file << generator;
	}
}

void TestCountSort()
{
	srand(time(0));
	std::vector<int32_t> array;
	for (int32_t i = 0; i < 1000; ++i)
		array.push_back(rand());

	std::vector<int32_t> out;
	
	std::cout << "The number before sorting are: " << std::endl;
	for (int32_t i = 0; i < 1000; ++i)
		std::cout << array[i] << " ";
	std::cout << std::endl << std::endl;
	CountingSort(array, out, SortOrder::Decreasing, 32768);
	std::cout << "The numbers after sorting are: " << std::endl;
	for (int32_t i = 0; i < 1000; ++i)
		std::cout << out[i] << " ";
	std::cout << std::endl << std::endl;
}
int main (int argc, char** argv)
{
	cout << std::endl;

	//VectorTest();
	//MatrixTest();
	//TestDCTs();
	//TestQuaternion();
	//TestHeap();
	//TestQuickSort();
	//Generate1MillionNumbers();
	//srand(time(0));
	//getchar();
	//TestQuickSort1();
	//TestCountingSort();
	//UseFindingSecondSmallestValue();
	//TestBST();
	//WatchStructureOfBST();
	//TestRedBlackTree();
	//WatchStructureOfRBT();
	//WatchStructureOfRBT();
	//WatchRBStructureOfRBT();
	//TestRedBlackTree();
	//TestInsertSorting();
	//TestAssignmentTime();
	//TestRandom();
	TestCountSort();
    return 0;
}