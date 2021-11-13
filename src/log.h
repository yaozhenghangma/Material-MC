#ifndef OUTLOG
#define OUTLOG

#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>

#include "MC_structure.h"

int WriteLog(Supercell &, MonteCarlo &, std::shared_ptr<spdlog::logger>);

#endif