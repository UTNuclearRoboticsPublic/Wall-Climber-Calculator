import yaml

def load_yaml(file_path):
	try:
		with open(file_path, 'r') as file:
			data = yaml.safe_load(file)
		return data
	except FileNotFoundError:
		print("File does not exist!")
		raise
	except:
		print("Load yaml broke!")
		raise 