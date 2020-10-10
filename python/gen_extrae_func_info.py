import os

func_names = [
	"spec_advance",
	"send_spec",
	"send_current",
	"wait_save_update_current",
	# "current_smooth",
	"yee_b",
	"yee_e",
	"get_part_seg_direction",
	"correct_coords",
	"send_emf_gc",
	# "emf_move_window",
	"wait_save_particles",
	"wait_save_emf_gc"
]

zpic_path = "../em2d-gaspi/zpic"

lines = []

for func_name in func_names:
	lines.append( os.popen(f"nm -a {zpic_path} | grep {func_name}").read().split("\n")[0].replace(" T ", "#") + "\n" ) 


with open("functions.dat", 'w') as file:
		file.writelines(lines)