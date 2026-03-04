#include "include/path_manager.h"

#include <sstream>
#include <stdexcept>

namespace easyparticle {

namespace {

// Normalize for stable comparisons/logs even when candidates contain "..".
std::filesystem::path NormalizePath(const std::filesystem::path &path) {
	return std::filesystem::weakly_canonical(path).lexically_normal();
}

// Resolve the running executable directory.
// On Linux we prefer /proc/self/exe, and fall back to current_path() if unavailable.
std::filesystem::path DetectExecDir() {
#if defined(__linux__)
	std::error_code ec;
	std::filesystem::path exe_file = std::filesystem::read_symlink("/proc/self/exe", ec);
	if (!ec && !exe_file.empty()) {
		return NormalizePath(exe_file).parent_path();
	}
#endif
	return std::filesystem::current_path();
}

} // namespace

// Cache executable directory once per process.
std::filesystem::path PathManager::ExecPath() {
	static const std::filesystem::path exec_path = DetectExecDir();
	return exec_path;
}

std::vector<std::filesystem::path> PathManager::CandidateDataDirs(
	const std::filesystem::path &exec_path
) {
	// Search order:
	// 1) Nearby runtime locations relative to executable.
	// 2) Standard share layouts for build/install trees.
	// 3) CMake-injected hints from parent projects.
	std::vector<std::filesystem::path> candidates = {
		exec_path,
		exec_path / "..",
		exec_path / "../share/easy-particle",
		exec_path / "../share/easyparticle",
		exec_path / "../../share/easy-particle",
		exec_path / "../../share/easyparticle",
		exec_path / "../../.."
	};

#ifdef EASY_PARTICLE_DATA_DIR
	candidates.emplace_back(EASY_PARTICLE_DATA_DIR);
#endif
#ifdef EASY_PARTICLE_SOURCE_DATA_DIR
	candidates.emplace_back(EASY_PARTICLE_SOURCE_DATA_DIR);
#endif

	return candidates;
}

std::filesystem::path PathManager::DataPath(const std::string &filename) {
	// Caller controls which data file to locate; PathManager only resolves directories.
	if (filename.empty()) {
		throw std::runtime_error("PathManager::DataPath() filename cannot be empty");
	}

	// Return the first directory that contains the requested file.
	const auto candidates = CandidateDataDirs(ExecPath());
	for (const auto &candidate : candidates) {
		const std::filesystem::path file_path = candidate / filename;
		if (std::filesystem::is_regular_file(file_path)) {
			return NormalizePath(candidate);
		}
	}

	std::ostringstream oss;
	oss << "Failed to locate data file '" << filename << "'. Checked paths:";
	for (const auto &candidate : candidates) {
		oss << "\n  - " << (candidate / filename);
	}
	throw std::runtime_error(oss.str());
}

} // namespace easyparticle
