#include <algorithm>
#include <iostream>
#include <sstream>
#include "pugixml/src/pugixml.hpp"

class XmlConfigurator
{
private:
	std::string _filename;
	int _startLevel;
	int _finishLevel;
	bool _isCpuComputeEnabled;
	bool _isOpenMPComputeEnabled;
	bool _isGpuComputeEnabled;
	bool _isCheckCorrectnessEnabled;
	bool _isNeedPrintInfo;

	void SetParams(const pugi::xml_document &);
	void SetDefaultStatus();

	std::string GetYesOrNo(bool condition)
	{
		if (condition)
			return "Yes";
		else
			return "No";
	}

	bool IntToBool(const int i)
	{
		return i != 0;
	}

public:

	XmlConfigurator()
	{
		_filename = "params.xml";
		SetDefaultStatus();
	}

	XmlConfigurator(std::string filename)
	{
		_filename = filename;
		SetDefaultStatus();
	}

	void LoadConfiguration();

	int GetStartLevel()
	{
		return _startLevel;
	}

	int GetFinishLevel()
	{
		return _finishLevel;
	}

	bool IsCpuComputeEnabled()
	{
		return _isCpuComputeEnabled;
	}

	bool IsOpenMPComputeEnabled()
	{
		return _isOpenMPComputeEnabled;
	}

	bool IsGpuComputeEnabled()
	{
		return _isGpuComputeEnabled;
	}

	bool IsCheckCorrectnessEnabled()
	{
		return _isCheckCorrectnessEnabled;
	}

	bool IsNeedPrintInfo()
	{
		return _isNeedPrintInfo;
	}

	void PrintConfiguration();
};