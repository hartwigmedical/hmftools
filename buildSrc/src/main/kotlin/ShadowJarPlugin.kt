import org.gradle.api.Plugin
import org.gradle.api.Project
import org.gradle.api.provider.Property
import org.gradle.api.Action
import org.gradle.kotlin.dsl.*

interface ShadowJarPluginExtension {
    val mainClass: Property<String>
}

// our own wrapper Gradle Shadow (https://github.com/johnrengelman/shadow) customise fat jar
class ShadowJarPlugin : Plugin<Project> {

    override fun apply(project: Project) {
        // apply the shawdow plugin
        project.plugins.apply("com.github.johnrengelman.shadow")

        // Add the 'buildJar' extension object
        val extension = project.extensions.create("shadowJar", ShadowJarPluginExtension::class.java)

        // get the archive version
        val versionPropertyName = "${project.name}.version"
        val projectVersion =
            if (project.hasProperty(versionPropertyName)) {
                project.property(versionPropertyName).toString()
            }
            else {
                project.logger.debug("Property $versionPropertyName missing. Using ${project.version}")
                project.version.toString()
            }

        // we need to use this such that project has a chance to configure before
        // we modify the shadow plugin
        // see https://stackoverflow.com/questions/57762850/how-to-use-gradle-plugin-configuration-in-the-plugin-apply
        project.afterEvaluate {
            project.tasks.withType<com.github.jengelman.gradle.plugins.shadow.tasks.ShadowJar> {
                logger.debug("main class = ${extension.mainClass.get()}")
                archiveVersion.set(projectVersion)
                manifest {
                    attributes["Main-Class"] = extension.mainClass.get()
                }
                archiveFileName.set("${project.name}-${project.version}-jar-with-dependencies.jar")

                // we want the build task to depend on this
                project.tasks.getByName("build").dependsOn(this)
            }
        }
    }
}

// https://stackoverflow.com/questions/62266477/gradle-custom-function-in-block-plugins
// this allows us to apply and configure the ShadowJar plugin like this:
//
// shadowJar {
//    mainClass.set("com.hartwig.hmftools.amber.AmberApplication")
// }
fun Project.shadowJar(configure: Action<ShadowJarPluginExtension>)
{
    this.apply<ShadowJarPlugin>()
    (this as org.gradle.api.plugins.ExtensionAware).extensions.configure("shadowJar", configure)
}
