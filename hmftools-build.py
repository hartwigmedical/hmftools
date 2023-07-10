import re
import subprocess
import sys
import xml.etree.ElementTree as ET

### Classes

class Maven:
    def __init__(self, pom_path):
        self.pom_path = pom_path
    
    def set_property(self, property, value):
        subprocess.run(['mvn', '-f', self.pom_path, 'versions:set-property',
                     '-DgenerateBackupPoms=false', f'-Dproperty={property}', f'-DnewVersion={value}'], check=True)
    
    def set_version(self, version):
        subprocess.run(['mvn', '-f', self.pom_path, 'versions:set',
                        '-DgenerateBackupPoms=false', f'-DnewVersion={version}'], check=True)
        
    def install(self, non_recursive=False):
        if non_recursive:
            subprocess.run(['mvn', '-N', '-f', self.pom_path, 'clean', 'install'], check=True)
        else:
            subprocess.run(['mvn', '-f', self.pom_path, 'clean', 'install'], check=True)

    def deploy(self, non_recursive=False):
        if non_recursive:
            subprocess.run(['mvn', '-N', '-f', self.pom_path, 'deploy', '-B'], check=True)
        else:
            subprocess.run(['mvn', '-f', self.pom_path, 'deploy', '-B'], check=True)


### Helper methods

def extract_hmftools_dependencies(pom_path):
    parsed_module_pom = ET.parse(pom_path)
    namespace = {'ns': 'http://maven.apache.org/POM/4.0.0'}
    dependencies = parsed_module_pom.root().findall('.//ns:dependencies/ns:dependency', namespace)
    hmftools_dependencies = set()
    for dep in dependencies:
        group_id = dep.find('ns:groupId', namespace).text
        artifact_id = dep.find('ns:artifactId', namespace).text
        if group_id == "com.hartwig":
            hmftools_dependencies.add(artifact_id)
    return hmftools_dependencies


### Script

# Check if a semver version is included as argument.
if len(sys.argv) != 2:
    print(f'Invalid arguments. Usage: {sys.argv[0]} <semver-version>')
    quit()

tag = sys.argv[1]

# check if the tag name is according to the regex
semver_pattern = '^([a-z-]+)-v?([0-9]+\.[0-9]+\.[0-9]+(?:-(?:alpha|beta)\.[0-9]+)?)$'
match = re.match(semver_pattern, tag)

if not match:
    print(f'Invalid tag (it does not match the regex pattern): {tag}')
    quit()

module = match.group(1)
version = match.group(2)

# parse all the hmftools modules the project depends on from the pom.xml
hmftools_dependencies = extract_hmftools_dependencies(f'{module}/pom.xml')

parent_pom = Maven('pom.xml')
module_pom = Maven(f'{module}/pom.xml')
dependencies_pom = [Maven(f'{hmf_dep}/pom.xml') for hmf_dep in hmftools_dependencies]

# Set versions in appropriate poms
# For the module we are targetting, we will use only the version part of the semver tag
# For all dependencies, we will use the entire semver tag 
parent_pom.set_property(f'{module}.version', version)
for hmf_dep in hmftools_dependencies:
    parent_pom.set_property(f'{hmf_dep}.version', tag)
parent_pom.set_version(tag)
module_pom.set_version(version)

# Build all submodules to check for build errors
parent_pom.install(non_recursive=True)
for dependency_pom in dependencies_pom:
    dependency_pom.install()
module_pom.install()

# Deploy if no errors
parent_pom.deploy(non_recursive=True)
for dependency_pom in dependencies_pom:
    dependency_pom.deploy()
module_pom.deploy()
