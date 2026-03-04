#ifndef __PATH_MANAGER_H__
#define __PATH_MANAGER_H__

#include <filesystem>
#include <string>
#include <vector>

namespace easyparticle {

class PathManager {
private:
	PathManager() = delete;
	~PathManager() = delete;
	PathManager(const PathManager&) = delete;
	PathManager& operator=(const PathManager&) = delete;

public:
	/// @brief Find directory containing the given data file
	/// @param filename data file name (e.g. amdc_ion_2020.txt)
	/// @returns directory containing filename
	static std::filesystem::path DataPath(const std::string &filename);

private:
	static std::filesystem::path ExecPath();
	static std::vector<std::filesystem::path> CandidateDataDirs(const std::filesystem::path &exec_path);
};

} // namespace easyparticle

#endif // __PATH_MANAGER_H__
