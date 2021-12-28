import org.gradle.kotlin.dsl.*
import java.time.Instant
import java.time.ZoneOffset
import java.time.format.DateTimeFormatter

/**
 * Generate the ${project.name}.version file in the resource directory
 */
plugins {
    java
}

tasks {
    val buildVersion by registering(Copy::class) {
        description = "Generate <project.name>.version resource file"

        val versionTemplateFile = File("${sourceSets.main.get().resources.srcDirs.first()}/${project.name}.version")

        if (!versionTemplateFile.exists())
        {
            throw GradleException("template file $versionTemplateFile missing")
        }

        from(versionTemplateFile)
        destinationDir = sourceSets.main.get().output.resourcesDir!!

        val props = mapOf("project.version" to project.version,
                                                "timestamp" to DateTimeFormatter
                                                    .ofPattern("yyyy-MM-dd HH:mm")
                                                    .withZone(ZoneOffset.UTC)
                                                    .format(Instant.now()))

        filter { line: String ->
            val regex = Regex("\\\$\\{([\\w.]+)}")
            regex.replace(line) { m ->
                val key = m.groupValues[1]
                if (key in props) props[key].toString() else m.value
            }
        }

        // always regenerate this file
        outputs.upToDateWhen { false }
    }

    classes {
        dependsOn(buildVersion)
    }

    // tell the normal resource processing to skip this file
    processResources {
        exclude("${project.name}.version")
    }
}
