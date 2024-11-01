import toml

# Load the existing pyproject.toml
pyproject_path = 'pyproject.toml'
pyproject = toml.load(pyproject_path)

# Read requirements from requirements.txt
with open('./src/POSST/requirements.txt') as f:
    requirements = f.read().splitlines()

# Update dependencies in the pyproject.toml
pyproject['project']['dependencies'] = requirements

# Write the updated pyproject.toml back to the file
with open(pyproject_path, 'w') as f:
    toml.dump(pyproject, f)

print("Updated pyproject.toml with dependencies from requirements.txt.")

