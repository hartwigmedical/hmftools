import re
import subprocess
import sys
import xml.etree.ElementTree as ET


def extract_hmftools_dependencies(dependencies):
    hmftools_dependencies = set()
    for dep in dependencies:
        group_id = dep.find('ns:groupId', namespace).text
        artifact_id = dep.find('ns:artifactId', namespace).text
        if group_id == "com.hartwig":
            hmftools_dependencies.add(artifact_id)
    return hmftools_dependencies





# Check if a semver version is included as argument.
if len(sys.argv) != 2:
    print(f'Invalid arguments. Usage: {sys.argv[0]} semver-version')
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
parsed_module_pom = ET.parse(f'{module}/pom.xml')
root = parsed_module_pom.getroot()
namespace = {'ns': 'http://maven.apache.org/POM/4.0.0'}
dependencies = root.findall('.//ns:dependencies/ns:dependency', namespace)
hmftools_dependencies = extract_hmftools_dependencies(dependencies)


# Set version of module in parent pom.
subprocess.run(
    ['mvn', '-f', 'pom.xml', 'versions:set-property', '-DgenerateBackupPoms=false', f'-Dproperty={module}.version',
     f'-DnewVersion={version}'],
    check=True)

# Set version of hmftools dependencies in parent pom.xml
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
