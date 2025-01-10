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
from typing import BinaryIO

import requests
import sys
import time
import jwt

from xml.etree import ElementTree
from argparse import ArgumentParser

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(stream=sys.stdout,
                    format='%(asctime)s [%(levelname)5s] - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)
logger.setLevel(logging.DEBUG)

SEMVER_REGEX = re.compile(
    r'^([a-z-]+)-v([0-9]+\.[0-9]+(?:\.[0-9]+)?(?:-(?:alpha|beta|rc)\.[0-9]+)?(?:_(?:[0-9a-zA-Z-]+(\.[0-9a-zA-Z-]+)*))?)$')

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
        logger.info('Starting Maven deploy_all')
        module_str = ','.join([m.name for m in modules])
        subprocess.run(['mvn', 'deploy', '-B', '-pl', module_str, '-am', '-DdeployAtEnd=true'], check=True)


class Docker:
    def __init__(self, module, version, is_external_release):
        self.module = module
        self.version = version
        self.is_external_release = is_external_release
        self.internal_image = f'europe-west4-docker.pkg.dev/hmf-build/hmf-docker/{self.module}:{self.version}'
        self.external_image = f'hartwigmedicalfoundation/{self.module}:{self.version}'

    def build(self):
        logger.info('Starting Docker build')
        with open("/workspace/entrypoint_template.sh", "r") as template, \
           open(f"/workspace/{self.module}/target/entrypoint.sh", "w") as output:
               for line in template:
                   output.write(line.rstrip().replace("__JAR_PATH__", f"/usr/share/java/{self.module}_v{self.version}.jar") + '\n')

        with open("/workspace/docker.sh", "w") as output:
            output.write('set -e\n')
            output.write(f'[ ! -f {self.module}/Dockerfile ] && echo "No Dockerfile for {self.module}" && exit 0\n')
            output.write('docker buildx create --name builder --driver docker-container --driver-opt network=cloudbuild --use\n')
            output.write(f'docker buildx build --add-host metadata.google.internal:169.254.169.254 {self.module} --load -t {self.internal_image} -t {self.external_image} --build-arg VERSION={self.version}\n')
            output.write(f'docker push {self.internal_image}\n')
            if self.is_external_release:
                output.write(f'cat /workspace/dockerhub.password | docker login -u hartwigmedicalfoundation --password-stdin\n')
                output.write(f'docker push {self.external_image}\n')
            else:
                output.write(f'echo Not pushing to external Docker registry as {self.version} is an internal-only version')


class GithubRelease:
    def __init__(self, tag_name: str, module: str, version: str, artifact_file: BinaryIO, private_key: str, github_client_id: str,
            github_installation_id: str):
        self.tag_name = tag_name
        self.module = module
        self.version = version
        self.artifact_file = artifact_file
        self.private_key = private_key
        self.github_client_id = github_client_id
        self.github_installation_id = github_installation_id
        self.release_name = f"{module} v{version}"

    def create(self):
        jwt = self._construct_jwt()
        token = self._refresh_token(jwt)
        logger.info("Successfully refreshed token")
        id = self._create_release(token)
        logger.info(f"Created release with id {id}")
        self._upload_artifacts(id, token)

    def _create_release(self, token: str):
        logger.info(f"Creating release [{self.release_name}]")
        request = {"tag_name": self.tag_name,
                "target_commitish": "master",
                "name": self.release_name,
                "body": f"Description of release {self.release_name}",
                "prerelease": True,
                "generate_release_notes": True
        }

        response = requests.post(self._construct_url("api"),
                json = request,
                headers = self._headers(token))
        logger.info(f"Response: {response.text}")
        response.raise_for_status()
        return response.json()["id"]

    def _upload_artifacts(self, id: str, token: str):
        logger.info(f"Uploading artifact to release {id}")
        headers = self._headers(token)
        headers["Content-Type"] = "application/octet-stream"
        base_url="{}/{}/assets?name".format(self._construct_url("uploads"), id)
        response = requests.post(f"{base_url}={self.module}_v{self.version}.jar", 
                headers = headers, 
                data = self.artifact_file.read())
        response.raise_for_status()
        logger.info(f"Uploaded {self.artifact_file.name}")

    def _construct_url(self, prefix):
        return f"https://{prefix}.github.com/repos/hartwigmedical/hmftools/releases"

    # This method gleaned from Github docs/examples
    def _construct_jwt(self):
        payload = {
            'iat': int(time.time()),
            # JWT expiration time (10 minutes maximum)
            'exp': int(time.time()) + 600,
            'iss': self.github_client_id
        }
        return jwt.encode(payload, self.private_key, algorithm='RS256')

    def _refresh_token(self, jwt: str):
        response = requests.post(f"https://api.github.com/app/installations/{self.github_installation_id}/access_tokens",
                headers = self._headers(jwt))
        response.raise_for_status()
        return response.json()["token"]

    def _headers(self, token: str):
        return {"Accept": "application/vnd.github+json",
                "Authorization": f"Bearer {token}",
                "X-GitHub-Api-Version": "2022-11-28"
        }


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
    parser.add_argument('github_key_path', help="Path to private key for the Github deployment bot")
    parser.add_argument('github_client_id', help="Client id for the deployment bot")
    parser.add_argument('github_installation_id', help="Installation id of the deployment bot")
    args = parser.parse_args()

    build_and_release(args.tag, args.github_key_path, args.github_client_id, args.github_installation_id)


def build_and_release(raw_tag: str, github_key: str, github_client_id: str, github_installation_id: str):
    logger.info('Starting build and release')
    match = SEMVER_REGEX.match(raw_tag)
    if not match:
        logger.info(f"Invalid tag: '{raw_tag}' (it does not match the regex pattern: '{SEMVER_REGEX.pattern}')")
        exit(1)
    module = match.group(1)
    version = match.group(2)
    is_external_release = re.compile(r'^[0-9\.]+(-rc\.[0-9]+)*$').match(version) != None

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

    Docker(module, version, is_external_release).build()
    if is_external_release:
        GithubRelease(raw_tag, module, version, open(f'/workspace/{module}/target/{module}-{version}-jar-with-dependencies.jar', 'rb'), 
            open(github_key, "r").read(), github_client_id, github_installation_id).create()
    else:
        logger.info(f'Skipping Github release creation as {version} is an internal-only version')

    logger.info('Complete build and release')

if __name__ == '__main__':
    main()
