#include "SceneManager.hpp"
#include "Core/AssetLoader.hpp"
#include "Core/Parser/OGEX.hpp"
#include "Core/Parser/DAE.hpp"

namespace Panda
{
	SceneManager::~SceneManager()
	{

	}

	int SceneManager::Initialize()
	{
		int result = 0;

		m_pScene = std::make_shared<Scene>();
		return result;
	}

	void SceneManager::Finalize()
	{

	}

	void SceneManager::Tick()
	{
		if (m_DirtyFlag)
			m_DirtyFlag = !(m_IsRenderingQueued);
	}

	int SceneManager::LoadScene(const char* sceneFileName)
	{
		// Judge file suffix.
		std::string str(sceneFileName);
		std::size_t found = str.find_last_of('.');
		std::string suffix = str.substr(found + 1);
		if (suffix == "ogex")
		{
			if (LoadOgexScene(sceneFileName))
			{
				m_pScene->LoadResource();
				m_DirtyFlag = true;
				m_IsRenderingQueued = false;
				return 0;
			}
			else
				return -1;
		}
		else if (suffix == "dae")
		{
			if (LoadDaeScene(sceneFileName))
			{
				m_pScene->LoadResource();
				m_DirtyFlag = true;
				m_IsRenderingQueued = false;
				return 0;
			}
			else
				return -1;
		}
		else
			return -1;
	}

	void SceneManager::ResetScene()
	{
		m_DirtyFlag = true;
	}

	bool SceneManager::LoadOgexScene(const char* ogexSceneFileName)
	{
		std::string ogexText = g_pAssetLoader->SyncOpenAndReadFileToString(ogexSceneFileName);

		if (ogexText.empty())
		{
			return false;
		}
		OgexParser ogexParser;
		m_pScene = ogexParser.Parse(ogexText);

		if (!m_pScene)
		{
			return false;
		}

		return true;
	}

	bool SceneManager::LoadDaeScene(const char* daeSceneFileName)
	{
		DaeParser daeParser;
		m_pScene = daeParser.Parse(daeSceneFileName);

		if (!m_pScene)
		{
			return false;
		}

		return true;
	}

	const Scene& SceneManager::GetScene()
	{
		return *m_pScene;
	}

	const Scene& SceneManager::GetSceneForPhysicalSimulation()
	{
		// TODO: we should perform CPU scene crop at here
		return *m_pScene;
	}

	bool SceneManager::IsSceneChanged()
	{
		return m_DirtyFlag;
	}

	void SceneManager::NotifySceneIsRenderingQueued()
	{
		m_IsRenderingQueued = true;
	}

	std::weak_ptr<SceneGeometryNode> SceneManager::GetSceneGeometryNode(std::string name)
	{
		auto it = m_pScene->LUTNameGeometryNode.find(name);
		if (it != m_pScene->LUTNameGeometryNode.end())
			return it->second;
		else 
			return std::weak_ptr<SceneGeometryNode>();
	}

	std::weak_ptr<SceneObjectGeometry> SceneManager::GetSceneGeometryObject(std::string key)
	{
		return m_pScene->Geometries.find(key)->second;
	}

	std::weak_ptr<BaseSceneNode> SceneManager::GetRootNode()
	{
		return m_pScene->SceneGraph;
	}
}
