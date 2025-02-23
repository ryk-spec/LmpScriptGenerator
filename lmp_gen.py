import os

class LAMMPSInputBuilder:
    # クラス変数：平衡化スケジュールの通し番号
    equil_counter = 0

    def __init__(self, do_minimization=True, data_prefix="hoge"):
        # コマンドリストを初期化
        self.commands = []
        self.max_skipped_index = 0
        # 初期設定
        self.initialize_coeff(data_prefix)
        self.data_prefix = data_prefix
        self.initialize_neighbors()
        self.initialize_thermo()
        
        if do_minimization:
            self.initialize_minimization()
    
    def initialize_coeff(self, data_prefix):
        """Class2力場用のpair_styleとスタイル（bond, angle, dihedral, improper）の設定
        （dataファイルにbond、angle、dihedral、improperなどのcoeffが記載されているため、pair_coeffは不要）
        """
        self.commands.append("units           real")
        self.commands.append("dimension       3")
        self.commands.append("boundary        p p p")
        self.commands.append("newton          on")
        self.commands.append("atom_style      full")
        self.commands.append("bond_style      class2")
        self.commands.append("angle_style     class2")
        self.commands.append("dihedral_style  class2")
        self.commands.append("improper_style  class2")
        self.commands.append("special_bonds   lj/coul 0 0 1")
        self.commands.append("kspace_style    pppm 1.0e-4")
        self.commands.append("pair_style      lj/class2/coul/long 12.0")
        self.commands.append("pair_modify     mix sixthpower")

        # read_dataコマンドを追加
        # hoge.${i}.dataが見つかった場合、hoge.${i+1}.restartを探索し、見つかったらそちらを読み込みます。
        # さらに、スキップする番号の最大値（i）をインスタンス変数に保持します。
        import glob, re
        files = glob.glob(f"{data_prefix}.*.data")
        max_num = -1
        selected_file = None
        pattern = r"^" + re.escape(data_prefix) + r"(?:\.(\d+))?\.data$"
        for f in files:
            match = re.match(pattern, f)
            if match:
                num = int(match.group(1)) if match.group(1) is not None else 0
                if num > max_num:
                    max_num = num
                    selected_file = f
        
        # スキップする番号の最大値をインスタンス変数に格納
        self.max_skipped_index = max_num
        
        # hoge.${i}.dataが見つかった場合、その番号i+1のrestartファイルを探索する
        if max_num >= 0:
            restart_candidate = f"{data_prefix}.{max_num+1}.restart"
            if os.path.exists(restart_candidate):
                selected_file = restart_candidate
        
        if selected_file is None:
            selected_file = f"{data_prefix}.data"
        if os.path.splitext(selected_file)[1] == ".data":
            self.commands.append(f"read_data {selected_file}")
        elif os.path.splitext(selected_file)[1] == ".restart":
            self.commands.append(f"read_restart {selected_file}")
        else:
            raise ValueError(f"Invalid data file: {selected_file}")
        
        if not os.path.exists(selected_file):
            raise ValueError(f"Data file not found: {selected_file}")
    
    def initialize_neighbors(self):
        """ペアリストの設定"""
        self.commands.append("neighbor 2.0 bin")
        self.commands.append("neigh_modify delay 0 every 1 check yes")
    
    def initialize_thermo(self, thermo_interval=1000):
        """
        thermo 出力の設定
        出力項目: ステップ、温度、ポテンシャルエネルギー、運動エネルギー、全エネルギー、圧力、ボリューム、密度
        密度は、compute と variable を用いて算出（dens = 総質量 / ボリューム）
        """
        self.commands.append("compute myMass all reduce sum mass")
        self.commands.append("variable dens equal c_myMass/vol")
        self.commands.append(f"thermo {thermo_interval}")
        self.commands.append("thermo_style custom step temp pe ke etotal press vol v_dens")
    
    def initialize_minimization(self):
        """
        エネルギー最小化の設定
        共役勾配法（cg）を用いた最小化
        """
        if LAMMPSInputBuilder.equil_counter <= self.max_skipped_index:
            return
        LAMMPSInputBuilder.equil_counter += 1
        self.commands.append("min_style cg")
        self.commands.append("min_modify dmax 0.1 line quadratic")
        self.commands.append("minimize 1.0e-4 1.0e-6 1000 10000")
        self.commands.append(f"write_data {self.data_prefix}.{LAMMPSInputBuilder.equil_counter}.data")
        self.commands.append(f"write_restart {self.data_prefix}.{LAMMPSInputBuilder.equil_counter}.restart")
    
    def add_nvt_equilibration(self, temp_start, temp_end, damping, run_steps,
                               timestep=1.0,
                               output_lmptrj=False,
                               trj_interval=1000):
        """
        nvt 平衡化スケジュールの追加
        :param temp_start: 初期温度
        :param temp_end: 最終温度
        :param damping: 温度制御用ダンピングパラメータ
        :param run_steps: 平衡化のステップ数
        :param timestep: 運動方程式の時間刻み
        :param output_lmptrj: （オプション）True の場合、lmptrj出力を行う
        :param trj_interval: （オプション）lmptrj出力間隔
        """
        LAMMPSInputBuilder.equil_counter += 1
        eq_id = LAMMPSInputBuilder.equil_counter

        # 運動方程式の時間刻みの設定
        self.commands.append(f"timestep {timestep}")

        # nvt 固有の fix を追加（通し番号付き）
        self.commands.append(f"fix {eq_id} all nvt temp {temp_start} {temp_end} {damping}")

        # lmptrj 出力が要求されている場合、軌道ダンプの設定を追加
        if output_lmptrj:
            self.commands.append(f"dump trj{eq_id} all custom {trj_interval} {self.data_prefix}.{eq_id}.lmptrj id type x y z")
            self.commands.append(f"dump_modify trj{eq_id} sort id")

        # restart 出力の設定とシミュレーションの分割実行
        self.commands.append("variable steps_per_ns equal 1000000")
        self.commands.append(f"variable rem equal {run_steps}")
        self.commands.append("variable ns equal 0")
        self.commands.append("label loop_start")
        self.commands.append("if \"${rem} >= ${steps_per_ns}\" then \"variable run_segment equal ${steps_per_ns}\"")
        self.commands.append("if \"${rem} < ${steps_per_ns}\" then \"variable run_segment equal ${rem}\"")
        self.commands.append("run ${run_segment}")
        self.commands.append(f"write_restart {self.data_prefix}.{eq_id}.restart")
        self.commands.append("variable ns equal ${ns}+1")
        self.commands.append("variable rem equal ${rem}-${run_segment}")
        self.commands.append("if \"${rem} > 0\" then \"jump SELF loop_start\"")

        # lmptrj 出力が要求されていた場合、ダンプの停止
        if output_lmptrj:
            self.commands.append(f"undump trj{eq_id}")

        # 平衡化後のデータ出力
        self.commands.append(f"write_data {self.data_prefix}.{eq_id}.data")
    
    def add_npt_equilibration(self, temp_start, temp_end, damping, pressure_start, pressure_end, pdamp, run_steps,
                               data_file=None):
        """
        npt 平衡化スケジュールの追加（nvt の形式に合わせたバージョン）
        :param temp_start: 初期温度
        :param temp_end: 最終温度
        :param damping: 温度制御用ダンピングパラメータ
        :param pressure_start: 初期圧力
        :param pressure_end: 最終温度
        :param pdamp: 圧力制御用ダンピングパラメータ
        :param run_steps: 平衡化のステップ数
        :param data_file: （オプション）出力する data file のパス。指定がない場合は自動生成します。
        """
        LAMMPSInputBuilder.equil_counter += 1
        eq_id = LAMMPSInputBuilder.equil_counter

        if data_file is None:
            data_file = f"{self.data_prefix}.{eq_id}.data"

        # 運動方程式の時間刻みの設定（npt でもデフォルトで 1.0 を使用）
        self.commands.append("timestep 1.0")
        # npt 固有の fix を追加（通し番号付き）
        self.commands.append(f"fix {eq_id} all npt temp {temp_start} {temp_end} {damping} iso {pressure_start} {pressure_end} {pdamp}")

        # シミュレーションの分割実行の設定
        self.commands.append("variable steps_per_ns equal 1000000")
        self.commands.append(f"variable rem equal {run_steps}")
        self.commands.append("variable ns equal 0")
        self.commands.append("label loop_start")
        self.commands.append("if \"${rem} >= ${steps_per_ns}\" then \"variable run_segment equal ${steps_per_ns}\"")
        self.commands.append("if \"${rem} < ${steps_per_ns}\" then \"variable run_segment equal ${rem}\"")
        self.commands.append("run ${run_segment}")
        self.commands.append(f"write_restart {self.data_prefix}.{eq_id}.restart")
        self.commands.append("variable ns equal ${ns}+1")
        self.commands.append("variable rem equal ${rem}-${run_segment}")
        self.commands.append("if \"${rem} > 0\" then \"jump SELF loop_start\"")

        # 平衡化後のデータ出力
        self.commands.append(f"write_data {data_file}")
    
    def add_command(self, command):
        """任意の LAMMPS コマンドを追加するメソッド"""
        self.commands.append(command)
    
    def generate_input(self):
        """コマンドリストからインプットスクリプトを生成"""
        return "\n".join(self.commands)
    
    def write_to_file(self, filename):
        """生成したインプットスクリプトをファイルに書き出す"""
        with open(filename, 'w') as f:
            f.write(self.generate_input())
# 使用例
if __name__ == "__main__":
    builder = LAMMPSInputBuilder(do_minimization=True, data_prefix="hoge")
    
    # nvt 平衡化スケジュール（data_file が存在する場合はスキップ、data_file は自動生成）
    builder.add_nvt_equilibration(
        temp_start=100.0, temp_end=300.0, damping=100.0, run_steps=5000
    )
    
    # npt 平衡化スケジュール（data_file が存在する場合はスキップ、data_file は自動生成）
    builder.add_npt_equilibration(
        temp_start=300.0, temp_end=300.0, damping=100.0,
        pressure_start=0.0, pressure_end=0.0, pdamp=1000.0, run_steps=5000
    )
        
    # 作成したインプットスクリプトの内容を表示
    print(builder.generate_input())
    
    # ファイルに書き出し
    builder.write_to_file("lammps_input.in")
