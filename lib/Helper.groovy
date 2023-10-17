class Helper {

    static def check_mandatory_parameter(Map params, String parameter_name){
        if ( !params[parameter_name]){
            println "You must specifiy a " + parameter_name
            System.exit(1)
        } else {
            return params[parameter_name]
        }
    }

    static def check_optional_parameters(Map params, List parameter_names){
        if (parameter_names.collect{element -> params[element]}.every{element -> element == false}){
            println "You must specifiy at least one of these options: " + parameter_names.join(", ")
            System.exit(1)
        }
    }

    static def check_parameter_value(String parameter_name, String value, List value_options){
        if (value_options.any{ it == value }){
            return value
        } else {
            println "The value for " + parameter_name + " must be one of " + value_options.join(", ")
            System.exit(1)
        }
    }

    static def complete_message(Map params, nextflow.script.WorkflowMetadata workflow, String version){
        // Display complete message
        println ""
        println "Ran the workflow: ${workflow.scriptName} ${version}"
        println "Command line    : ${workflow.commandLine}"
        println "Completed at    : ${workflow.complete}"
        println "Duration        : ${workflow.duration}"
        println "Success         : ${workflow.success}"
        println "Exit status     : ${workflow.exitStatus}"
        println "Work directory  : ${workflow.workDir}"
        println ""
        println "Parameters used:"
        println "----------------"
        params.each{ k, v ->
            if (v){
                println "${k}: ${v}"
            }
        }
        println "====================================================="
    }

    static def success_message() {
        println ""
        println "====================================================="
        println "            WORKFLOW SUCCESSFULLY COMPLETED"
        println "====================================================="
    }

    static def error_message(nextflow.script.WorkflowMetadata workflow){
        // Display error message
        println ""
        println "====================================================="
        println "            WORKFLOW FAILED"
        println "====================================================="
        println "Workflow execution stopped with the following message:"
        println "  " + workflow.errorMessage
    }

}
