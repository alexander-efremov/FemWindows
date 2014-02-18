enum ModelDataType
{ Model11,
	Model21,
	Model41,
	Model81,
	Model161,
	Model321,
	Model641,
	Model1281
};

class ModelDataProvider
{
public:

	double* GetModelData(enum ModelDataType type);
};