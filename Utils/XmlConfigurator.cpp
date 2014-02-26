#include "XmlConfigurator.h"

void XmlConfigurator::LoadConfiguration()
{
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(_filename.c_str());
	SetParams(doc);
}

void XmlConfigurator::SetParams(const pugi::xml_document& doc)
{
	
	pugi::xml_node startup = doc.child("Startup");
	_startLevel = atoi(startup.child("StartLevel").child_value());
	_finishLevel = atoi(startup.child("FinishLevel").child_value());
	_isCpuComputeEnabled = IntToBool(atoi(startup.child("IsCpuComputeEnabled").child_value()));
	_isOpenMPComputeEnabled = IntToBool(atoi(startup.child("IsOpenMPComputeEnabled").child_value()));
	_isGpuComputeEnabled = IntToBool(atoi(startup.child("IsGpuComputeEnabled").child_value()));
	_isCheckCorrectnessEnabled = IntToBool(atoi(startup.child("IsCheckCorrectnessEnabled").child_value()));
	_isNeedPrintInfo = IntToBool(atoi(startup.child("IsPrintingEnabled").child_value()));
}

void XmlConfigurator::PrintConfiguration()
{
	std::cout << std::endl;
	std::cout << "Start level: " << _startLevel << std::endl;
	std::cout << "Finish level: " << _finishLevel << std::endl;
	std::cout << "Is CPU compute enabled: " << GetYesOrNo(_isCpuComputeEnabled) << std::endl;
	std::cout << "Is OpenMP compute enabled: " << GetYesOrNo(_isOpenMPComputeEnabled) << std::endl;
	std::cout << "Is GPU compute enabled: " << GetYesOrNo(_isGpuComputeEnabled) << std::endl;
	std::cout << "Is error checking enabled: " << GetYesOrNo(_isCheckCorrectnessEnabled) << std::endl;
	std::cout << "Is printing enabled: " << GetYesOrNo(_isNeedPrintInfo) << std::endl;
	std::cout << std::endl;
}

void XmlConfigurator::SetDefaultStatus()
{
	_startLevel = 0;
	_finishLevel = 0;
	_isCpuComputeEnabled = false;
	_isOpenMPComputeEnabled = false;
	_isGpuComputeEnabled = false;
	_isNeedPrintInfo = false;
}