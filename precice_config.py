
### !! NOT UP TO DATE !! ###

def header():
    text = "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n" \
        + "<precice-configuration>\n"\
        + "  <log>\n"\
        + "    <sink\n"\
        + "      filter=\"%Severity% > debug\"\n"\
        + "      format=\"---[precice] %ColorizedSeverity% %Message%\"\n"\
        + "      enabled=\"true\" />\n"\
        + "  </log>\n\n"\
        + "  <solver-interface dimensions=\"2\">\n"\
        
    return text

def add_cell(n, no_cells, dt, t_end):
    n_str = str(n)
    text =  "    <data:scalar name=\"cell_internal_heat_" + n_str + "\" />\n\n"\
        + "    <mesh name=\"electrostatics_mesh_" + n_str + "\">\n"\
        + "      <use-data name=\"cell_internal_heat_" + n_str + "\" />\n"\
        + "    </mesh>\n\n"\
        + "    <participant name=\"electrostatics_" + n_str + "\">\n"\
        + "      <use-mesh name=\"electrostatics_mesh_" + n_str + "\" provide=\"yes\" />\n"\
        + "      <write-data name=\"cell_internal_heat_" + n_str + "\" mesh=\"electrostatics_mesh_" + n_str + "\" />\n"\
        + "    <mesh name=\"thermal_mesh_" + n_str + "\">\n"\
        + "      <use-data name=\"cell_internal_heat_" + n_str + "\" />\n"\
        + "    </mesh>\n\n"
    
    #if n == (no_cells - 1):
    #
    #elif n == 0:
    #
    #else: 
    #    text +=  "    </participant>\n\n"\
    #        + "    <participant name=\"thermal_" + n_str + "\">\n"\
    #        + "      <use-mesh name=\"electrostatics_mesh_" + n_str + "\" from=\"electrostatics_" + n_str + "\" />\n"\
    #        + "      <use-mesh name=\"thermal_mesh_" + n_str + "\" provide=\"yes\" />\n"\
    #        + "      <mapping:nearest-neighbor\n"\
    #        + "        direction=\"read\"\n"\
    #        + "        from=\"electrostatics_mesh_" + n_str + "\"\n"\
    #        + "        to=\"thermal_mesh_" + n_str + "\"\n"\
    #        + "        constraint=\"consistent\" />\n"\
    #        + "      <read-data name=\"cell_internal_heat_" + n_str + "\" mesh=\"thermal_mesh_" + n_str + "\" /> \n"\
    #        + "    </participant>\n\n"

    
    
    text += "    <m2n:sockets from=\"electrostatics_" + n_str + "\" to=\"thermal_" + n_str + "\" exchange-directory=\"..\" />\n\n"\
        + "    <coupling-scheme:serial-explicit>\n"\
        + "      <participants first=\"electrostatics_" + n_str + "\" second=\"thermal_" + n_str + "\" />\n"\
        + "      <time-window-size value=\"" + str(dt) + "\" />\n"\
        + "      <max-time value=\"" + str(t_end) + "\" />\n"\
        + "      <exchange data=\"cell_internal_heat_" + n_str + "\" mesh=\"electrostatics_mesh_" + n_str + "\" from=\"electrostatics_" + n_str + "\" to=\"thermal_" + n_str + "\" />\n"\
        + "    </coupling-scheme:serial-explicit>\n"\
        
    return text


def closing():
    text = "  </solver-interface>\n"\
        + "</precice-configuration>\n"
    
    return text



def adapter_config(no_cells):
    for i in range(no_cells):
        i_str = str(i)
        text = "{\n" \
            + "  \"participant_name\": \"thermal_" + i_str + "\",\n" \
            + "  \"config_file_name\": \"../../precice-config.xml\",\n" \
            + "  \"interface\": {\n" \
            + "      \"coupling_mesh_name\": \"thermal_mesh_" + i_str + "\",\n" \
            + "      \"write_data_name\": \"cell_to_cell_heat_" + i_str + "\",\n" \
            + "      \"read_data_name\": \"cell_to_cell_temp_" + i_str + "\"\n" \
            + "    }\n" \
            + "}\n"
        
        filename = "thermal/precice_adapter_configs/prenics_adapter_cell_" + i_str + ".json"

        with open(filename, "w") as file:
            file.writelines(text)
            
    

def build_config(no_cells, dt, t_end):

    text = header()
    
    for i in range(no_cells):
        text += add_cell(i, no_cells, dt, t_end)

    text += closing()


    with open("precice-config.xml", "w") as file:
        file.writelines(text)



