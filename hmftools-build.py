"""
This script does the following:
- Whenever a user pushes a git tag in the format '<tool>-<version>' this script will parse that tag.
  The 'v' character prefixing the version will be stripped.
- The script will build the tool with the given version and deploy it using Maven with version equal to the 'version'
  part of the tag.
- The script will build all direct hmf-tools dependencies of the given tool (i.e. 'hmf-common') and deploy those with
  the version set equal to '<tool>-<version>'.

By building and deploying each tool together with their dependencies, we ensure each deployment is isolated.
Example:
    - 'git tag neo-v1.0.0' will start building the 'neo' tool and deploy it with version '1.0.0'.
    - 'neo' has a dependency on 'hmf-common', so it will deploy 'hmf-common' with version 'neo-1.0.0'.
    - The pom of the deployed 'neo' will be updated such that the 'hmf-common' dependency will point towards the correct
      version.
"""
import re
import subprocess
from xml.etree import ElementTree
from argparse import ArgumentParser

SEMVER_REGEX = re.compile(
    r'^([a-z-]+)-v([0-9]+\.[0-9]+(?:\.[0-9]+)?(?:-(?:alpha|beta)\.[0-9]+)?(?:_(?:[0-9a-zA-Z-]+(\.[0-9a-zA-Z-]+)*))?)$')


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


def extract_hmftools_dependencies(pom_path):
    namespace = {'ns': 'http://maven.apache.org/POM/4.0.0'}
    # First, obtain a list of all modules defined in the parent
    parsed_parent_pom = ElementTree.parse('pom.xml')
    modules = parsed_parent_pom.getroot().findall('.//ns:modules/ns:module', namespace)
    module_set = {module.text for module in modules}
    # Then, obtain dependencies on these modules from target module
    parsed_module_pom = ElementTree.parse(pom_path)
    dependencies = parsed_module_pom.getroot().findall('.//ns:dependencies/ns:dependency', namespace)
    hmftools_dependencies = set()
    for dep in dependencies:
        group_id = dep.find('ns:groupId', namespace).text
        artifact_id = dep.find('ns:artifactId', namespace).text
        if group_id == "com.hartwig" and artifact_id in module_set:
            hmftools_dependencies.add(artifact_id)
    return hmftools_dependencies


def main():
    parser = ArgumentParser(
        description="A tool for automatically building and deploying individual modules in HMF-tools.")
    parser.add_argument('tag', help="The semantic versioning tag in the following format: <tool-name>-<version>")
    args = parser.parse_args()

    build_and_release(args.tag)


def build_and_release(raw_tag: str):
    match = SEMVER_REGEX.match(raw_tag)
    if not match:
        print(f"Invalid tag: '{raw_tag}' (it does not match the regex pattern: '{SEMVER_REGEX.pattern}')")
        exit(1)
    module = match.group(1)
    version = match.group(2)
    # Clean the raw_tag such that it only includes the groups captured by the regex
    # For example: raw_tag = orange-v1.0.0 then tag = orange-1.1.0
    tag = f'{module}-{version}'

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


if __name__ == '__main__':
    main()
