#include <iostream>
#include <map>
#include <string>

class BMI_Var
{
private:
	std::string name;
	std::string type;
	std::string units;
	bool set;
	bool get;
public:
	// methods
	BMI_Var()
	{
		this->set = false;
		this->get = false;
	}
	BMI_Var(std::string name_in, std::string type_in, std::string units_in,
		bool set_in, bool get_in)
	{
		this->name = name_in;
		this->type = type_in;
		this->units = units_in;
		this->set = set_in;
		this->get = get_in;
	};

	std::string GetName() { return this->name; }
	std::string GetType() { return this->type; };
	std::string GetUnits() { return this->units; }
	void SetUnits(std::string units_in) { this->units = units_in; };
	bool GetSet(void) { return this->set; };
	void SetSet(bool in_in) { this->set = in_in; };
	bool GetGet(void) { return this->get; };
	void SetGet(bool out_in) { this->get = out_in; };
};
