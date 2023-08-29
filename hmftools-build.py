import re
import subprocess
import sys
import xml.etree.ElementTree as ET

### Classes

class Maven:
    def __init__(self, pom_path, name=''):
        self.pom_path = pom_path
        self.name = name
    
    def set_property(self, property, value):
        subprocess.run(['mvn', '-f', self.pom_path, 'versions:set-property',
                     '-DgenerateBackupPoms=false', f'-Dproperty={property}', f'-DnewVersion={value}'], check=True)
    
    def set_version(self, version):
        subprocess.run(['mvn', '-f', self.pom_path, 'versions:set',
                        '-DgenerateBackupPoms=false', f'-DnewVersion={version}'], check=True)

    @staticmethod
    def deploy_all(*modules):
        module_str = ','.join([m.name for m in modules])
        subprocess.run(['mvn', 'deploy', '-B', '-pl', module_str, '-am', '-DdeployAtEnd=true'])


### Helper methods

def extract_hmftools_dependencies(pom_path):
    namespace = {'ns': 'http://maven.apache.org/POM/4.0.0'}

    # First, obtain a list of all modules defined in the parent
    parsed_parent_pom = ET.parse('pom.xml')
    modules = parsed_parent_pom.getroot().findall('.//ns:modules/ns:module', namespace)
    module_set = {module.text for module in modules}
    # Then, obtain dependencies on these modules from target module
    parsed_module_pom = ET.parse(pom_path)
    dependencies = parsed_module_pom.getroot().findall('.//ns:dependencies/ns:dependency', namespace)
    hmftools_dependencies = set()
    for dep in dependencies:
        group_id = dep.find('ns:groupId', namespace).text
        artifact_id = dep.find('ns:artifactId', namespace).text
        if group_id == "com.hartwig" and artifact_id in module_set:
            hmftools_dependencies.add(artifact_id)
    return hmftools_dependencies


## Script

# Check if a semver version is included as argument.
if len(sys.argv) != 2:
    print(f'Invalid arguments. Usage: {sys.argv[0]} <semver-version>')
    quit()

tag = sys.argv[1]

# check if the tag name is according to the regex
semver_pattern = '^([a-z-]+)-(v?[0-9]+\.[0-9]+(?:\.[0-9]+)?(?:-(?:alpha|beta)\.[0-9]+)?(?:_(?:[0-9a-zA-Z-]+(\.[0-9a-zA-Z-]+)*))?)$'
match = re.match(semver_pattern, tag)

if not match:
    print(f'Invalid tag (it does not match the regex pattern): {tag}')
    quit()

module = match.group(1)
version = match.group(2)

# parse all the hmftools modules the project depends on from the pom.xml
hmftools_dependencies = extract_hmftools_dependencies(f'{module}/pom.xml')

parent_pom = Maven('pom.xml')
module_pom = Maven(f'{module}/pom.xml', name=module)
dependencies_pom = [Maven(f'{hmf_dep}/pom.xml', name=hmf_dep) for hmf_dep in hmftools_dependencies]

# Set versions in appropriate poms
# For the module we are targeting, we will use only the version part of the semver tag
# For all dependencies, we will use the entire semver tag
parent_pom.set_property(f'{module}.version', version)
for hmf_dep in hmftools_dependencies:
    parent_pom.set_property(f'{hmf_dep}.version', tag)
parent_pom.set_version(tag)
module_pom.set_version(version)

Maven.deploy_all(module_pom, *dependencies_pom)
