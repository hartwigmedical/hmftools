import org.gradle.api.Plugin
import org.gradle.api.Project
import org.gradle.api.Action
import org.gradle.api.tasks.SourceSetContainer
import org.gradle.api.tasks.bundling.Jar
import org.gradle.kotlin.dsl.*

interface TestsJarPluginExtension

// build a test jar for the test classes
// in our code we have some test classes inside the test directory that are used by other projects
// see https://stackoverflow.com/questions/5644011/multi-project-test-dependencies-with-gradle/60138176#60138176
//
class TestsJarPlugin : Plugin<Project> {

    override fun apply(project: Project) {
        // Add the extension object
        project.extensions.create("testsJar", TestsJarPluginExtension::class.java)

        // build a test jar for the test classes
        // in our code we have some test classes inside the test directory that are used by other projects
        // see https://stackoverflow.com/questions/5644011/multi-project-test-dependencies-with-gradle/60138176#60138176
        val testsJar by project.tasks.registering(Jar::class) {
            archiveClassifier.set("tests")
            from(project.extensions.getByType(SourceSetContainer::class)["test"].output)
        }

        project.configurations.create("tests")

        project.artifacts.add("tests", testsJar)
    }
}

// https://stackoverflow.com/questions/62266477/gradle-custom-function-in-block-plugins
// this allows us to apply and configure the TestsJar plugin like this:
//
// testsJar {
// }
fun Project.testsJar(configure: Action<TestsJarPluginExtension>)
{
    this.apply<TestsJarPlugin>()
    (this as org.gradle.api.plugins.ExtensionAware).extensions.configure("testsJar", configure)
}