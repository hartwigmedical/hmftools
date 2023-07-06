import sys
import re
import subprocess
import xml.etree.ElementTree as ET

if len(sys.argv) != 2:
    print(f'Invalid arguments. Usage: {sys.argv[0]} semver-version')
    quit()

tag = sys.argv[1]
# check if the tag name is according to the regex
semver_pattern = '^(\w+)-v?((\d+)(\.\d+){0,2})(-(alpha|beta)\.[0-9]+)?$'
match = re.match(semver_pattern, tag)

if not match:
    print(f'Invalid tag (it does not match the regex pattern): {tag}')
    quit()

module = match.group(1)


# parse all the pom modules the project depends on from the pom.xml
parsed_pom = ET.parse(f'{module}/pom.xml')
root = parsed_pom.getroot()

namespace = {'ns': 'http://maven.apache.org/POM/4.0.0'}

dependencies = root.findall('.//ns:dependencies/ns:dependency', namespace)

hmftools_dependencies = set()

for dep in dependencies:
    # Extract relevant information from each dependency
    group_id = dep.find('ns:groupId', namespace).text
    artifact_id = dep.find('ns:artifactId', namespace).text
    if group_id == "com.hartwig":
        hmftools_dependencies.add(artifact_id)



# Set version of module in parent pom.
subprocess.run(['mvn', '-f', 'pom.xml', 'versions:set-property', '-DgenerateBackupPoms=false', f'-Dproperty={module}.version',  f'-DnewVersion={tag}'],
               check=True)

# Set version of dependencies in parent pom.xml
for hmf_dep in hmftools_dependencies:
    subprocess.run(['mvn', '-f', 'pom.xml', 'versions:set-property', '-DgenerateBackupPoms=false',
                    f'-Dproperty={hmf_dep}.version', f'-DnewVersion={tag}'], check=True)

# Set version of parent pom
subprocess.run(['mvn', '-f', 'pom.xml', 'versions:set',
               '-DgenerateBackupPoms=false', f'-DnewVersion={tag}'], check=True)

# Build all submodules to check for build errors
subprocess.run(['mvn', '-N', '-f', '"pom.xml"',
               'clean', 'install'], check=True)
for hmf_dep in hmftools_dependencies:
    subprocess.run(
        ['mvn', '-f', f'{hmf_dep}/pom.xml', 'clean', 'install'], check=True)
subprocess.run(
    ['mvn', '-f', f'{module}/pom.xml', 'clean', 'install'], check=True)

# Deploy if no errors
subprocess.run(['mvn', '-N', '-f', 'pom.xml', 'deploy', '-B'], check=True)
for hmf_dep in hmftools_dependencies:
    subprocess.run(
        ['mvn', '-f', f'{hmf_dep}/pom.xml', 'deploy', '-B'], check=True)

subprocess.run(['mvn', '-f', f'{module}/pom.xml', 'deploy', '-B'], check=True)
